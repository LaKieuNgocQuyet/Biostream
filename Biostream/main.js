const { app, BrowserWindow, ipcMain } = require('electron/main');
const path = require('node:path');
const fs = require('fs');

let win;

function createWindow () {
  win = new BrowserWindow({
    frame: false,
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      nodeIntegration: false,
      sandbox: false  
    }
  });

  win.setMenu(null);
  win.loadFile('index.html');
  win.webContents.openDevTools();
  win.maximize();
}

// Xử lý các sự kiện điều khiển cửa sổ
ipcMain.on('minimize-window', () => {
  if (win) win.minimize();
});

ipcMain.on('maximize-window', () => {
  if (win) {
    if (win.isMaximized()) {
      win.unmaximize();
    } else {
      win.maximize();
    }
  }
});

ipcMain.on('close-window', () => {
  if (win) win.close();
});

app.whenReady().then(() => {
  createWindow();

  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) {
      createWindow();
    }
  });
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});
