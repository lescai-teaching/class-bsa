const grid = document.querySelector("#material-grid");
const count = document.querySelector("#catalog-count");
const viewerTitle = document.querySelector("#viewer-title");
const viewerFrame = document.querySelector("#viewer-frame");
const openDirect = document.querySelector("#open-direct");

const state = {
  materials: [],
  selected: null
};

function cardTypeLabel(item) {
  if (item.type === "html") return "Pagina HTML";
  if (item.type === "pdf") return "Documento PDF";
  return "Contenuto non rilevato";
}

function setSelected(item) {
  state.selected = item;

  document.querySelectorAll(".material-card").forEach((card) => {
    card.toggleAttribute("aria-current", card.dataset.name === item.name);
  });

  viewerTitle.textContent = item.title || item.name;
  viewerFrame.className = "viewer-frame";
  viewerFrame.innerHTML = "";

  if (!item.href) {
    openDirect.hidden = true;
    viewerFrame.classList.add("empty");
    viewerFrame.innerHTML = "<p>La cartella non contiene un file index.html o un PDF visualizzabile.</p>";
    return;
  }

  openDirect.href = item.href;
  openDirect.hidden = false;

  const frame = document.createElement("iframe");
  frame.title = item.type === "pdf" ? `PDF: ${item.title || item.name}` : `Pagina: ${item.title || item.name}`;
  frame.src = item.href;
  viewerFrame.append(frame);
}

function createCard(item) {
  const card = document.createElement("button");
  card.type = "button";
  card.className = "material-card";
  card.dataset.name = item.name;

  const title = document.createElement("span");
  title.className = "card-title";
  title.textContent = item.title || item.name;

  const meta = document.createElement("span");
  meta.className = "card-meta";
  meta.textContent = cardTypeLabel(item);

  if (item.file) {
    const file = document.createElement("span");
    file.className = "card-file";
    file.textContent = item.file;
    card.append(title, meta, file);
  } else {
    card.append(title, meta);
  }

  card.addEventListener("click", () => {
    setSelected(item);
    history.replaceState(null, "", `#${encodeURIComponent(item.name)}`);
  });

  return card;
}

function renderMaterials(materials) {
  grid.innerHTML = "";

  if (materials.length === 0) {
    count.textContent = "Nessuna cartella pubblicata.";
    grid.innerHTML = '<div class="empty-card">Aggiungi cartelle con un file index.html o PDF a questo ramo.</div>';
    return;
  }

  const foldersLabel = materials.length === 1 ? "cartella pubblicata" : "cartelle pubblicate";
  count.textContent = `${materials.length} ${foldersLabel}`;
  materials.forEach((item) => grid.append(createCard(item)));

  const hashName = decodeURIComponent(window.location.hash.replace(/^#/, ""));
  const requested = materials.find((item) => item.name === hashName);
  setSelected(requested || materials[0]);
}

async function loadMaterials() {
  try {
    const response = await fetch("materials.json", { cache: "no-store" });
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    const data = await response.json();
    state.materials = Array.isArray(data.items) ? data.items : [];
    renderMaterials(state.materials);
  } catch (error) {
    count.textContent = "Indice non disponibile.";
    grid.innerHTML = '<div class="empty-card">Impossibile caricare materials.json.</div>';
    viewerFrame.classList.add("empty");
    viewerFrame.innerHTML = "<p>Rigenera l'indice o verifica la pubblicazione della pagina.</p>";
  }
}

loadMaterials();
