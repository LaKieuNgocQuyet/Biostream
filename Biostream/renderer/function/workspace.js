

document.addEventListener("DOMContentLoaded", () => {
	const tabBar = document.querySelector(".tab-bar");
	const tabContent = document.querySelector(".tab-content");
	const analysisPanel = document.querySelector(".analysis-panel");
	const analysisButtons = document.querySelectorAll(".analysis-panel button");
	const tabs = document.querySelectorAll(".tab-bar .tab");

	tabs.forEach((tab, index) => {
		tab.dataset.order = index;
	});

	function createNewTab(icon, name) {
		const tabName = `${name}`;

		// Hủy active
		document.querySelectorAll(".tab").forEach(tab => tab.classList.remove("active"));
		document.querySelectorAll(".tab-pane").forEach(pane => pane.classList.remove("active"));

		// Tạo tab
		const newTab = document.createElement("div");
		const newIcon = document.createElement("i");
		newIcon.className = icon.className;
		newTab.appendChild(newIcon);

		newTab.className = "tab active";


		const textNode = document.createTextNode(" " + tabName);
		newTab.appendChild(textNode);
		newTab.addEventListener("click", () => activateTab(newTab));
		tabBar.appendChild(newTab);

		// Tạo nội dung tab
		const newContent = document.createElement("div");
		newContent.className = "tab-pane active";
		fetch('./renderer/components/variant_calling_tab/parameters.html')
			.then(response => {
				if (!response.ok) {
				throw new Error('Network response was not ok');
				}
				return response.text();
			})
			.then(html => {
				newContent.innerHTML = html;
			})
			.catch(error => {
				console.error('Lỗi:', error);
  		});

		tabContent.appendChild(newContent);
	}

	function activateTab(tabElement) {
		document.querySelectorAll(".tab").forEach(tab => tab.classList.remove("active"));
		document.querySelectorAll(".tab-pane").forEach(pane => pane.classList.remove("active"));

		tabElement.classList.add("active");

		const index = [...tabBar.children].indexOf(tabElement);
		const panes = document.querySelectorAll(".tab-pane");
		if (panes[index]) {
			panes[index].classList.add("active");
		}
	}

	analysisButtons.forEach(button => {
		button.addEventListener("click", () => {
			const name = button.textContent.trim();
			const analysisBtn = document.getElementById("analysis-button");
			const icon = button.querySelector("i");
			createNewTab(icon, name);
			analysisPanel.classList.remove("show");
			analysisBtn.classList.remove("onclick");
		});
	});

});
