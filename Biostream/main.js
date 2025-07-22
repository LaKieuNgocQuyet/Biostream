

const openPanel = document.getElementById('openPanel');
const closePanel = document.getElementById('closePanel');
const sidePanel = document.getElementById('sidePanel');
const overlay = document.getElementById('overlay');

openPanel.addEventListener('click', () => {
  sidePanel.classList.add('show');
  overlay.classList.remove('hidden');
});

closePanel.addEventListener('click', () => {
  sidePanel.classList.remove('show');
  overlay.classList.add('hidden');
});

overlay.addEventListener('click', () => {
  sidePanel.classList.remove('show');
  overlay.classList.add('hidden');
});
