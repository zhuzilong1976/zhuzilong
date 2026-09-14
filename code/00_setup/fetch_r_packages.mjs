import fs from "node:fs/promises";
import path from "node:path";

const RVER = process.env.RVER || "4.6";
const repo = `https://cran.r-project.org/bin/windows/contrib/${RVER}`;
const outDir = "work/rpkgs";
await fs.mkdir(outDir, { recursive: true });

// 随 R 一起安装的基础包/推荐包，不需要下载
const BASE = new Set([
  "base","compiler","datasets","grDevices","graphics","grid","methods","parallel",
  "splines","stats","stats4","tcltk","tools","utils","MASS","Matrix","boot","class",
  "cluster","codetools","foreign","KernSmooth","lattice","mgcv","nlme","nnet","rpart",
  "spatial","survival","R"
]);

async function fetchText(url) {
  const r = await fetch(url, { headers: { "User-Agent": "codex-research" } });
  if (!r.ok) throw new Error(`${url} -> ${r.status}`);
  return r.text();
}

console.log(`CRAN binary repo: ${repo}`);
const packagesTxt = await fetchText(`${repo}/PACKAGES`);
await fs.writeFile(`${outDir}/PACKAGES`, packagesTxt);

// 解析 PACKAGES（Debian-control 风格，字段可能有续行）
const records = new Map();
let cur = null;
let lastKey = null;
for (const rawLine of packagesTxt.split(/\r?\n/)) {
  const line = rawLine.replace(/\r$/, "");
  if (line.trim() === "") {
    if (cur && cur.Package) records.set(cur.Package, cur);
    cur = null;
    lastKey = null;
    continue;
  }
  if (/^\s/.test(line)) {
    // 续行：DCF 格式，接到上一个字段
    if (cur && lastKey) cur[lastKey] = `${cur[lastKey]} ${line.trim()}`;
    continue;
  }
  const m = line.match(/^([A-Za-z0-9._-]+):\s?(.*)$/);
  if (m && cur !== null || (m && cur === null)) {
    cur = cur || {};
    cur[m[1]] = m[2];
    lastKey = m[1];
  }
}
if (cur && cur.Package) records.set(cur.Package, cur);
console.log(`PACKAGES 解析到 ${records.size} 个包`);

function depsOf(rec) {
  const out = new Set();
  for (const field of ["Depends", "Imports", "LinkingTo"]) {
    const v = rec[field];
    if (!v) continue;
    for (const part of v.split(",")) {
      const name = part.trim().split(/[\s(]/)[0];
      if (name && name !== "R" && !BASE.has(name)) out.add(name);
    }
  }
  return [...out];
}

const roots = (process.env.ROOTS || "scTenifoldKnk,scTenifoldNet,locfdr,fgsea,preprocessCore")
  .split(",")
  .map((s) => s.trim())
  .filter(Boolean);

const missingRoots = roots.filter((r) => !records.has(r));
const needed = new Map();
const queue = roots.filter((r) => records.has(r));
const unresolved = new Set();

while (queue.length) {
  const name = queue.shift();
  if (needed.has(name)) continue;
  const rec = records.get(name);
  if (!rec) {
    unresolved.add(name);
    continue;
  }
  needed.set(name, rec);
  for (const d of depsOf(rec)) {
    if (!needed.has(d)) queue.push(d);
  }
}

console.log(`\n根包: ${roots.join(", ")}`);
if (missingRoots.length) console.log(`CRAN 二进制仓库中没有: ${missingRoots.join(", ")}`);
if (unresolved.size) console.log(`无法解析的依赖: ${[...unresolved].join(", ")}`);
console.log(`需要安装的包共 ${needed.size} 个:`);
for (const [n, r] of [...needed].sort()) {
  console.log(`  ${n} ${r.Version} (compilation=${r.NeedsCompilation || "no"})`);
}

// 拓扑排序：被依赖者先装
const order = [];
const state = new Map();
function visit(name, stack = []) {
  if (state.get(name) === 2) return;
  if (state.get(name) === 1) return; // 循环依赖，忽略
  state.set(name, 1);
  const rec = needed.get(name);
  if (rec) for (const d of depsOf(rec)) if (needed.has(d)) visit(d, [...stack, name]);
  state.set(name, 2);
  order.push(name);
}
for (const n of needed.keys()) visit(n);

await fs.writeFile(`${outDir}/install_order.txt`, order.join("\n"));
console.log(`\n安装顺序已写入 work/rpkgs/install_order.txt（${order.length} 个）`);

// 下载
const failed = [];
let downloaded = 0;
let totalBytes = 0;
for (const name of order) {
  const rec = needed.get(name);
  const file = `${name}_${rec.Version}.zip`;
  const dest = path.join(outDir, file);
  try {
    await fs.access(dest);
    continue;
  } catch {}
  try {
    const r = await fetch(`${repo}/${file}`, { headers: { "User-Agent": "codex-research" } });
    if (!r.ok) throw new Error(`HTTP ${r.status}`);
    const buf = Buffer.from(await r.arrayBuffer());
    await fs.writeFile(dest, buf);
    downloaded++;
    totalBytes += buf.length;
    console.log(`  downloaded ${file} (${(buf.length / 1048576).toFixed(2)} MB)`);
  } catch (e) {
    failed.push(`${file}: ${e.message}`);
  }
}

console.log(`\n下载完成: ${downloaded} 个新文件, ${(totalBytes / 1048576).toFixed(1)} MB`);
if (failed.length) {
  console.log("下载失败:");
  for (const f of failed) console.log("  " + f);
}
