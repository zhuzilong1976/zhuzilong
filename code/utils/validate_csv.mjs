import fs from "node:fs";

function parse(text) {
  const rows = [];
  let row = [];
  let field = "";
  let inQuotes = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (inQuotes) {
      if (c === '"') {
        if (text[i + 1] === '"') {
          field += '"';
          i++;
        } else inQuotes = false;
      } else field += c;
    } else if (c === '"') {
      inQuotes = true;
    } else if (c === ",") {
      row.push(field);
      field = "";
    } else if (c === "\n") {
      row.push(field);
      rows.push(row);
      row = [];
      field = "";
    } else if (c !== "\r") {
      field += c;
    }
  }
  if (field !== "" || row.length) {
    row.push(field);
    rows.push(row);
  }
  return rows;
}

const files = [
  "outputs/templates/dataset_survey.csv",
  "outputs/templates/gene_panel.csv",
  "outputs/templates/metadata_schema.csv",
  "outputs/templates/disease_options.csv"
];

for (const file of files) {
  const raw = fs.readFileSync(file, "utf8");
  const text = raw.replace(/^\uFEFF/, "");
  const rows = parse(text).filter((r) => r.length > 1 || (r.length === 1 && r[0] !== ""));
  const n = rows[0].length;
  let bad = 0;
  rows.forEach((r, i) => {
    if (r.length !== n) {
      bad++;
      console.log(`  BAD ${file} line ${i + 1}: fields=${r.length} expected=${n} :: ${r.slice(0, 3).join(" | ")}`);
    }
  });
  console.log(`${file} -> dataRows=${rows.length - 1} cols=${n} badRows=${bad} header=${rows[0].slice(0, 3).join(",")}`);
  fs.writeFileSync(file, "\uFEFF" + text, "utf8");
}
console.log(`BOM restored on ${files.length} files.`);
