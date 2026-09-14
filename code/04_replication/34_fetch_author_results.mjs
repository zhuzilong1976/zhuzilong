const files = [
  ["作者 TREM2 主结果", "inst/manuscript/TREM2/Results/drTREM2.csv"],
  ["作者稳定性 r1", "inst/manuscript/STABILITY/r1sko_Trem2.csv"],
  ["作者方向性结果", "inst/manuscript/TREM2/dirResults/dirTREM2.csv"]
];

const base = "https://raw.githubusercontent.com/cailab-tamu/scTenifoldKnk/master/";
const clean = (s) => s.replace(/"/g, "");

for (const [name, path] of files) {
  const r = await fetch(base + path);
  const lines = (await r.text()).trim().split(/\r?\n/);
  const hdr = lines[0].split(",").map(clean);
  const iG = hdr.indexOf("gene");
  const iP = hdr.indexOf("p.adj");
  const iD = hdr.indexOf("distance");
  const iZ = hdr.indexOf("Z");
  console.log(`\n=== ${name} ===`);
  console.log(`列: ${hdr.join(" | ")}`);

  let sig = 0;
  const ds = [];
  for (let i = 1; i < lines.length; i++) {
    const f = lines[i].split(",").map(clean);
    if (f.length < 5) continue;
    if (parseFloat(f[iP]) < 0.05) sig++;
    ds.push([f[iG], parseFloat(f[iD]), parseFloat(f[iZ])]);
  }
  ds.sort((a, b) => b[1] - a[1]);
  const med = ds[Math.floor(ds.length / 2)][1];
  console.log(`基因数 = ${ds.length}   p.adj<0.05 = ${sig} 个 (${((100 * sig) / ds.length).toFixed(2)}%)`);
  console.log(`距离最大 5: ${ds.slice(0, 5).map((x) => `${x[0]}=${x[1].toExponential(2)}`).join(", ")}`);
  console.log(`中位距离 = ${med.toExponential(2)}   最大/中位 = ${(ds[0][1] / med).toFixed(1)}`);
}
