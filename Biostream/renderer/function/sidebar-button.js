document.addEventListener("DOMContentLoaded", () => {
  const sidebar = document.querySelector(".sidebar");
  const toggleBtn = document.querySelector(".toggle-btn");
  const contentContainer = document.querySelector(".content-container");
  const toggleIcon = toggleBtn.querySelector("i");
  const analysisBtn = document.getElementById("analysis-button");
  const analysisPanel = document.querySelector(".analysis-panel");
  const sidebarbutton = document.querySelector(".sidebar button");

  let isCollapsed = false;
  let hideTimeout = null;
  let isAnalysisPanelOpen = false; // 👉 Biến theo dõi trạng thái panel

  const showToggleButton = () => {
    if (isAnalysisPanelOpen) return; // 👉 Không hiện nếu panel đang mở
    clearTimeout(hideTimeout);
    toggleBtn.classList.add("show");
  };

  const hideToggleButton = () => {
    clearTimeout(hideTimeout);
    hideTimeout = setTimeout(() => {
      toggleBtn.classList.remove("show");
    }, 3000);
  };

  // Hover hiện/ẩn toggleBtn
  sidebar.addEventListener("mouseenter", showToggleButton);
  toggleBtn.addEventListener("mouseenter", showToggleButton);
  sidebar.addEventListener("mouseleave", hideToggleButton);
  toggleBtn.addEventListener("mouseleave", hideToggleButton);

  // Toggle sidebar
  toggleBtn.addEventListener("click", () => {
    isCollapsed = !isCollapsed;
    sidebar.classList.toggle("collapsed", isCollapsed);
    toggleBtn.classList.toggle("collapsed", isCollapsed);
    contentContainer.classList.toggle("collapsed", isCollapsed);
    toggleIcon.classList.toggle("bi-chevron-double-left", !isCollapsed);
    toggleIcon.classList.toggle("bi-chevron-double-right", isCollapsed);
    updateAnalysisPanelPosition();
  });

  // Click nút phân tích
  analysisBtn.addEventListener("click", (e) => {
    e.stopPropagation();
    analysisPanel.classList.toggle("show");
    toggleBtn.classList.remove("show");
    isAnalysisPanelOpen = analysisPanel.classList.contains("show"); // Cập nhật trạng thái
    updateAnalysisPanelPosition();
  });

  // Click ngoài panel để ẩn
  document.addEventListener("click", (e) => {
    if (!analysisPanel.contains(e.target) && e.target !== analysisBtn) {
      analysisPanel.classList.remove("show");
      analysisBtn.classList.remove("onclick");
      isAnalysisPanelOpen = false; // Reset lại
    }
  });

  // Hiệu ứng nhấn sidebar button
  sidebarbutton.addEventListener("click", () => {
    sidebarbutton.classList.toggle("onclick");
  });


  function updateAnalysisPanelPosition() {
  const sidebarWidth = sidebar.offsetWidth;
  analysisPanel.style.left = `${sidebarWidth}px`;
  }

});
