import fs from "node:fs/promises";

const base = "https://raw.githubusercontent.com/cailab-tamu/scTenifoldKnk/master/";
const files = [
  "inst/manuscript/AHR/Results/Preenterocytes.RData",
  "inst/manuscript/CFTR/Results/SRS3161261.RData",
  "inst/manuscript/CFTR/Results/SRS4245406.RData",
  "inst/manuscript/DMD/Results/GSM4116571.RData",
  "inst/manuscript/HNF4A-HNF4G/Results/GSM3477499.RData",
  "inst/manuscript/MALAT1/Results/MALAT1.RData",
  "inst/manuscript/MECP2/Results/SRS3059998.RData",
  "inst/manuscript/MECP2/Results/SRS3059999.RData",
  "inst/manuscript/NKX2-1/Results/GSM3716703.RData",
  "inst/manuscript/NEGCONTROL/SRS3161261_Pulmonary alveolar type II cellsAkap7.RData"
];

await fs.mkdir("work/author_networks", { recursive: true });
let ok = 0, bytes = 0;
for (const f of files) {
  const name = f.split("/").pop();
  const dest = `work/author_networks/${name}`;
  try {
    await fs.access(dest);
    console.log(`已存在 ${name}`);
    continue;
  } catch {}
  try {
    const r = await fetch(base + f, { redirect: "follow" });
    if (!r.ok) throw new Error(`HTTP ${r.status}`);
    const b = Buffer.from(await r.arrayBuffer());
    await fs.writeFile(dest, b);
    bytes += b.length;
    ok++;
    console.log(`${name}  ${(b.length / 1048576).toFixed(2)} MB`);
  } catch (e) {
    console.log(`FAIL ${name}: ${e.message}`);
  }
}
console.log(`\n完成：新下载 ${ok} 个，共 ${(bytes / 1048576).toFixed(1)} MB`);
