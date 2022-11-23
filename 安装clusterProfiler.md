## 安装 clusterProfiler 时遇到的一些问题

### systemfont 的报错

```bash
fatal error: fontconfig/fontconfig.h: No such file or directory compilation terminated.
```

安装

```bash
sudo apt install libfontconfig1-dev
```

### Rcrul 的报错

```bash
checking for curl-config... no
Cannot find curl-config
ERROR: configuration failed for package ‘RCurl’
```

安装

```bash
sudo apt install libcurl3-dev
```

### XML 的报错

```bash
Cannot find xml2-config
```

安装

```bash
sudo apt install libxml2-dev
```

### openssl 的报错

```bash
Configuration failed because openssl was not found.
```

安装

```bash
sudo apt install libssl-dev
```

