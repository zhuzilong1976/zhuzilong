const series = ["GSE279180", "GSE279183", "GSE301908", "GSE277430", "GSE279181"];

function ftpPath(gse) {
  const n = gse.match(/GSE(\d+)/)[1];
  const prefix = `GSE${n.slice(0, Math.max(1, n.length - 3))}nnn`;
  return `https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}/${gse}/suppl/`;
}

for (const gse of series) {
  const url = ftpPath(gse);
  console.log(`\n===== ${gse} =====`);
  console.log(url);
  try {
    const r = await fetch(url, { headers: { "User-Agent": "codex-research" } });
    if (!r.ok) {
      console.log(`  HTTP ${r.status}`);
      continue;
    }
    const html = await r.text();
    // 解析文件名与大小
    const rows = [...html.matchAll(/<a href="([^"]+)">[^<]*<\/a>\s*<\/td>\s*<td[^>]*>([^<]*)<\/td>\s*<td[^>]*>([^<]*)<\/td>/g)];
    if (rows.length === 0) {
      // 退回到简单匹配
      const files = [...html.matchAll(/href="([^"]+\.(?:gz|zip|tar|h5|h5ad|rds|RDS|mtx|tsv|csv))"/g)].map((m) => m[1]);
      for (const f of new Set(files)) console.log(`  ${f}`);
    } else {
      for (const m of rows) {
        const name = m[1];
        if (name.startsWith("/") || name.startsWith("..") || name.includes("?")) continue;
        console.log(`  ${name}  |  ${m[2].trim()}  |  ${m[3].trim()}`);
      }
    }
  } catch (e) {
    console.log(`  ERR ${e.message}`);
  }
}
