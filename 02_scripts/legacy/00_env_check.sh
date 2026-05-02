#!/bin/bash
echo "========== Pillar 2 环境检查 =========="
echo "时间: $(date)"
R -e "cat(\"R版本:\", R.version.string, \"\\n\"); library(Seurat); cat(\"Seurat:\", packageVersion(\"Seurat\"), \"\\n\")"
