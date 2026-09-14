// 用 Node 解析并下载 Python 依赖树（cp310 / win_amd64 / py3-none-any），供 pip 离线安装
import fs from "node:fs/promises";

const PY = "3.10";
const PLAT = "win_amd64";
const outDir = "work/wheels";
await fs.mkdir(outDir, { recursive: true });

const UA = { "User-Agent": "codex-research" };
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

// 不用安装的包（Windows 难以编译 / 本流程用不到）
const SKIP = new Set(["velocyto", "gimmemotifs", "jupyter", "jupyterlab",
                      "notebook"]);

// 显式版本偏好（CellOracle 的 pin）
const PREFERRED = {
  numpy: "==1.26.4",
  pandas: "<=1.5.3",
  matplotlib: "<3.7",
  anndata: "<=0.10.8",
  numba: ">=0.50.1,<0.61",
  scanpy: "<1.11"
};

const cache = new Map();
async function getJSON(url, tries = 3) {
  if (cache.has(url)) return cache.get(url);
  for (let i = 0; i < tries; i++) {
    try {
      const r = await fetch(url, { headers: UA });
      if (!r.ok) throw new Error(`HTTP ${r.status}`);
      const j = await r.json();
      cache.set(url, j);
      return j;
    } catch (e) {
      if (i === tries - 1) throw e;
      await sleep(1500);
    }
  }
}

function cmpVer(a, b) {
  const pa = a.split(/[.\-+]/).map((x) => (/^\d+$/.test(x) ? +x : 0));
  const pb = b.split(/[.\-+]/).map((x) => (/^\d+$/.test(x) ? +x : 0));
  for (let i = 0; i < Math.max(pa.length, pb.length); i++) {
    const d = (pa[i] || 0) - (pb[i] || 0);
    if (d) return d > 0 ? 1 : -1;
  }
  return 0;
}

function isPrerelease(v) {
  // 去掉 .postN 后，只要还含字母（a/b/rc/dev 等）即视为预发布
  return /[a-zA-Z]/.test(v.replace(/\.post\d+$/i, ""));
}

// 解析简单算子集合
function parseSpec(spec) {
  return spec
    .split(",")
    .map((s) => s.trim())
    .filter(Boolean)
    .map((s) => {
      const m = s.match(/^(===|==|>=|<=|~=|!=|<|>)\s*(.+)$/);
      return m ? { op: m[1], v: m[2].trim() } : null;
    })
    .filter(Boolean);
}

function satisfies(v, ops) {
  return ops.every(({ op, v: t }) => {
    const c = cmpVer(v, t);
    switch (op) {
      case "==": case "===": return c === 0;
      case "!=": return c !== 0;
      case ">=": return c >= 0;
      case "<=": return c <= 0;
      case ">": return c > 0;
      case "<": return c < 0;
      case "~=": {
        const parts = t.split(".").map(Number);
        const upper = parts.slice(0, -1);
        upper[upper.length - 1] += 1;
        return c >= 0 && cmpVer(v, upper.join(".")) < 0;
      }
      default: return true;
    }
  });
}

function wheelOK(filename) {
  if (!filename.endsWith(".whl")) return false;
  if (/(musllinux|manylinux|macosx|aarch64|arm64|i686|win32\.whl)/.test(filename)) return false;
  if (filename.includes("py3-none-any") || filename.includes("py2.py3-none-any")) return true;
  if (!filename.includes(PLAT)) return false;
  return filename.includes(`cp${PY.replace(".", "")}`) || filename.includes("abi3") ||
         filename.includes("py3-none");
}

// 依赖标记求值：针对 python 3.10 / win32
function markerTrue(marker) {
  if (!marker) return true;
  let s = marker;
  s = s.replace(/python_version/g, `"${PY}"`);
  s = s.replace(/python_full_version/g, `"${PY}.0"`);
  s = s.replace(/sys_platform/g, `"win32"`);
  s = s.replace(/platform_system/g, `"Windows"`);
  s = s.replace(/os_name/g, `"nt"`);
  s = s.replace(/platform_machine/g, `"AMD64"`);
  s = s.replace(/platform_python_implementation/g, `"CPython"`);
  s = s.replace(/implementation_name/g, `"cpython"`);
  s = s.replace(/extra/g, `""`);
  s = s.replace(/\band\b/g, "&&").replace(/\bor\b/g, "||").replace(/\bnot\b/g, "!");
  try {
    // 仅允许字面量与比较运算
    if (!/^[\s"'0-9.!=<>&|()a-zA-Z_\-+*]*$/.test(s)) return true;
    return Function(`"use strict";return (${s});`)() === true;
  } catch {
    return true;   // 无法解析时保守地纳入
  }
}

function parseReq(req) {
  const [main, marker] = req.split(";").map((x) => x.trim());
  const m = main.match(/^([A-Za-z0-9_.\-]+)\s*(?:\[[^\]]*\])?\s*(.*)$/);
  if (!m) return null;
  return { name: m[1].toLowerCase().replace(/_/g, "-"),
           spec: m[2] ? parseSpec(m[2]) : [], marker };
}

async function pickVersion(name, ops) {
  const j = await getJSON(`https://pypi.org/pypi/${name}/json`);
  const vers = Object.keys(j.releases).filter((v) => !isPrerelease(v) && j.releases[v].length);
  const cand = vers.filter((v) => satisfies(v, ops)).sort(cmpVer).reverse();
  for (const v of cand) {
    const files = j.releases[v].filter((f) => {
      if (!wheelOK(f.filename) || f.yanked) return false;
      // 关键：过滤掉不支持当前 Python 版本的 wheel
      if (f.requires_python) {
        try {
          if (!satisfies(PY, parseSpec(f.requires_python))) return false;
        } catch { return false; }
      }
      return true;
    });
    if (files.length) {
      const f = files.find((x) => x.filename.includes(PLAT)) || files[0];
      return { version: v, file: f };
    }
  }
  return null;
}

const roots = process.argv.slice(2);
const chosen = new Map();
const queue = [];

for (const r of roots) {
  const p = parseReq(r);
  if (p) queue.push(p);
}

while (queue.length) {
  const req = queue.shift();
  if (SKIP.has(req.name)) continue;
  const pref = PREFERRED[req.name];
  const ops = pref ? [...parseSpec(pref), ...req.spec] : req.spec;
  const cur = chosen.get(req.name);
  if (cur) {
    if (ops.length && !cur.ops.length) { cur.ops = ops; }
    else continue;
  }
  let pick;
  try {
    pick = await pickVersion(req.name, ops);
  } catch (e) {
    console.log(`  !! ${req.name}: ${e.message}`);
    continue;
  }
  if (!pick) { console.log(`  !! 找不到可用 wheel: ${req.name} ${JSON.stringify(ops)}`); continue; }
  chosen.set(req.name, { version: pick.version, file: pick.file, ops });
  console.log(`${req.name}==${pick.version}  (${(pick.file.size / 1048576).toFixed(1)} MB)`);

  // 取该版本的依赖
  let vj;
  try {
    vj = await getJSON(`https://pypi.org/pypi/${req.name}/${pick.version}/json`);
  } catch { continue; }
  for (const rd of vj.info.requires_dist || []) {
    const p = parseReq(rd);
    if (!p) continue;
    if (!markerTrue(p.marker)) continue;
    if (!chosen.has(p.name)) queue.push(p);
  }
  await sleep(60);
}

console.log(`\n解析完成：${chosen.size} 个包`);
let total = 0;
const manifest = [];
for (const [name, c] of chosen) {
  const dest = `${outDir}/${c.file.filename}`;
  try {
    await fs.access(dest);
  } catch {
    const r = await fetch(c.file.url, { headers: UA });
    const b = Buffer.from(await r.arrayBuffer());
    await fs.writeFile(dest, b);
    total += b.length;
    console.log(`  下载 ${c.file.filename} (${(b.length / 1048576).toFixed(1)} MB)`);
  }
  manifest.push({ name, version: c.version, file: c.file.filename });
}
await fs.writeFile(`${outDir}/manifest.json`, JSON.stringify(manifest, null, 2));
console.log(`新增下载 ${(total / 1048576).toFixed(1)} MB，清单已写入 work/wheels/manifest.json`);
