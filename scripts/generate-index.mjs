import { access, readdir, stat, writeFile } from "node:fs/promises";
import path from "node:path";

const root = process.cwd();
const supportDirectories = new Set([".git", ".github", "assets", "scripts", "node_modules"]);

function collatorSort(a, b) {
  return a.localeCompare(b, "it", { numeric: true, sensitivity: "base" });
}

function encodePath(...segments) {
  return segments.map((segment) => encodeURIComponent(segment)).join("/");
}

function titleFromFolder(folderName) {
  return folderName
    .split("_")
    .filter(Boolean)
    .map((word) => word.charAt(0).toUpperCase() + word.slice(1).toLowerCase())
    .join(" ");
}

async function exists(filePath) {
  try {
    await access(filePath);
    return true;
  } catch {
    return false;
  }
}

async function getFolderEntry(folderName) {
  const folderPath = path.join(root, folderName);
  const files = await readdir(folderPath, { withFileTypes: true });
  const indexPath = path.join(folderPath, "index.html");

  if (await exists(indexPath)) {
    return {
      name: folderName,
      title: titleFromFolder(folderName),
      type: "html",
      href: `${encodePath(folderName)}/`,
      file: "index.html"
    };
  }

  const pdf = files
    .filter((entry) => entry.isFile() && entry.name.toLowerCase().endsWith(".pdf"))
    .map((entry) => entry.name)
    .sort(collatorSort)[0];

  if (pdf) {
    return {
      name: folderName,
      title: titleFromFolder(folderName),
      type: "pdf",
      href: encodePath(folderName, pdf),
      file: pdf
    };
  }

  return {
    name: folderName,
    title: titleFromFolder(folderName),
    type: "unsupported",
    href: null,
    file: null
  };
}

async function main() {
  const entries = await readdir(root, { withFileTypes: true });
  const folders = [];

  for (const entry of entries) {
    if (!entry.isDirectory()) continue;
    if (supportDirectories.has(entry.name)) continue;
    if (entry.name.startsWith(".")) continue;

    const folderStats = await stat(path.join(root, entry.name));
    if (!folderStats.isDirectory()) continue;
    folders.push(entry.name);
  }

  const items = await Promise.all(folders.sort(collatorSort).map(getFolderEntry));
  const manifest = {
    generatedAt: new Date().toISOString(),
    items
  };

  await writeFile(
    path.join(root, "materials.json"),
    `${JSON.stringify(manifest, null, 2)}\n`,
    "utf8"
  );
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
