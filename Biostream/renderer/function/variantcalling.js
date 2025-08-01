const { exec } = require('child_process');

document.getElementById('variantcalling-btn').addEventListener('click', () => {
    exec('wsl bash /path/to/your_script.sh', (error, stdout, stderr) => {
        if (error) {
            console.error(`Error: ${error.message}`);
            return;
        }
        if (stderr) {
            console.error(`Stderr: ${stderr}`);
            return;
        }
        console.log(`✅ Output:\n${stdout}`);
    });
});
