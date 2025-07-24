const cpuContent = document.querySelector('.cpu-content');
const memoryContent = document.querySelector('.memory-content');
const diskContent = document.querySelector('.disk-info span');

async function updateMonitoringData() {
  const memory = window.electronAPI.getSystemMemoryInfo();
  const usedMemory = memory.total - memory.free;
  const memoryUsage = (usedMemory / memory.total) * 100;
  memoryContent.textContent = `${memoryUsage.toFixed(2)}%`;


  const cpu = window.electronAPI.getCPUUsage();
  const percent = cpu.percentCPUUsage * 100;
  cpuContent.textContent = `${percent.toFixed(2)}%`;



  const disks = await window.electronAPI.getDiskInfo();
  if (disks.length > 0) {
    const disk = disks[0]; // ổ đầu tiên
    const usedGB = (disk.used / 1e9).toFixed(2);
    const totalGB = (disk.total / 1e9).toFixed(2);
    diskContent.textContent = `Disk: ${usedGB} GB used (limit ${totalGB} GB)`;
  }
}

setInterval(updateMonitoringData, 700);
