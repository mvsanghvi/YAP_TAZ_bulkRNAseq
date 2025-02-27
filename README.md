# YAP_TAZ_bulkRNAseq
Bulk RNA sequencing of WT, YAP KO, and YAP/TAZ KO H9 cells

## Getting started

```bash
docker build -t yaptaz .
docker run -it --name ytcontainer yaptaz
```

To run STAR
```bash
./STAR --help
```
On 64 GB Codespace, takes 40 mins to make the genome index file