import fs from "node:fs/promises";

const libs = {
  KEGG_2021_Human: "KEGG_2021_Human",
  Reactome_2022: "Reactome_2022",
  GO_Biological_Process_2021: "GO_Biological_Process_2021",
  MSigDB_Hallmark_2020: "MSigDB_Hallmark_2020"
};

await fs.mkdir("work/gmt", { recursive: true });

for (const [name, lib] of Object.entries(libs)) {
  const url = `https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=${lib}`;
  try {
    const r = await fetch(url, { headers: { "User-Agent": "codex-research" } });
    if (!r.ok) throw new Error(`HTTP ${r.status}`);
    const text = await r.text();
    const nTerms = text.trim().split(/\r?\n/).length;
    await fs.writeFile(`work/gmt/${name}.gmt`, text);
    console.log(`${name}: ${nTerms} 个基因集, ${(text.length / 1048576).toFixed(2)} MB`);
  } catch (e) {
    console.log(`FAIL ${name}: ${e.message}`);
  }
}
