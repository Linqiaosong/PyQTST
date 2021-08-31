## Installation

You can install PyQTST by using ```pip``` or ```pip3```.

such as:

```bash
pip install PyQTST
```

or

```bash
pip3 install PyQTST
```

You can also install PyQTST from Github repo.

such as:

```bash
git clone https://github.com/Linqiaosong/PyQTST
cd PyQTST
python setup.py install
```

Using source code to install is also available.

such as:

```bash
tar -zxvf PyQTST-2.*.tar.gz
cd PyQTST-2.*
python setup.py install
```

## API for Shermo

**Shermo Website**: [http://sobereva.com/soft/shermo](http://sobereva.com/soft/shermo)

Shermo Citation: T. Lu, Q. Chen, Shermo: A general code for calculating molecular thermodynamic properties, ChemRxiv (2020), DOI: 10.26434/chemrxiv.12278801

We strongly recommend adding Shermo path to ```environment variables```.

In Windows System, you can change your ```environment variables``` at ```Control Panel\System and Security\System```

In Linux System, you can adding following codes to ```~/.bashrc```

```bash
export PATH=$PATH:/your/shermo/path
export Shermopath=/your/shermo/path
```