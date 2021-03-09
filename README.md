## Getting started

0. Make sure you have access to ACI:
   * Follow the account setup procedure [here](https://www.icds.psu.edu/computing-services/account-setup/)

1. Make sure you access to ```/gpfs/group/jtw13/default```
    * You can email [Jason Wright](mailto:astrowright@gmail.com) to be added

2. Set up an SSH key for GitHub (if you haven't already):
   * Follow the Linux directions [here](https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)

3. Clone this repo:
```bash
cd ~/work/
git clone git@github.com:Tusay/psu_seti_576.git
```

4. Activate the correct Anaconda environment:
    * For use on ACI:
    ```bash
    source /gpfs/group/jtw13/default/anaconda3/bin/activate
    conda activate seti_aci
    ```
    * For use on CyberLamp (e.g., with GPU)
    ```bash
    source /gpfs/group/jtw13/default/anaconda3/bin/activate
    conda activate seti_gpu
    ```

## Updating turbo_seti

```bash
cd /gpfs/group/jtw13/default/turbo_seti/
git pull
cd ..
pip install ./turbo_seti
```
