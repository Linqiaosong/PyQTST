# PyQTST

```
   d8888b. db    db  .d88b.  d888888b .d8888. d888888b 
   88  `8D `8b  d8' .8P  Y8. `~~88~~' 88'  YP `~~88~~' 
   88oodD'  `8bd8'  88    88    88    `8bo.      88    
   88~~~      88    88    88    88      `Y8b.    88 
   88         88    `8P  d8'    88    db   8D    88    
   88         YP     `Y88'Y8    YP    `8888Y'    YP    
```


## Latest Version 2.0.7

PyQTST is a software package for calculating chemical reaction rate constant by using transition state theory, which is written by Python3. It provides transition state theory calculation of mono-molecular and bi-molecular reactions, and also provides Skodje-Truhlar method to calculate one-dimensional quantum tunneling transmission coefficient that will be used in correcting rate constant. An API for Shermo is also provided by PyQTST.

**Github Website**: [https://github.com/Linqiaosong/PyQTST](https://github.com/Linqiaosong/PyQTST)

**PyPI Website**: [https://pypi.org/project/PyQTST](https://pypi.org/project/PyQTST)

**Online Documents**: [https://github.com/Linqiaosong/PyQTST/wiki](https://github.com/Linqiaosong/PyQTST/wiki)

**PyQTST Citation**: Q. Lin, PyQTST (version), https://github.com/Linqiaosong/PyQTST, (year)

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
