const { contextBridge, ipcRenderer} = require('electron');
const si = require('systeminformation');
const os = require('os');

contextBridge.exposeInMainWorld('electronAPI', {
	getSystemMemoryInfo: () => process.getSystemMemoryInfo(),
	getCPUUsage: () => process.getCPUUsage(),
	getDiskInfo: async () => {
    try {
      const data = await si.fsSize();
      return data.map(disk => ({
        mount: disk.mount,
        used: disk.used,
        total: disk.size,
      }));
    } catch (err) {
      console.error('Disk info error:', err);
      return [];
    }
  }
});
