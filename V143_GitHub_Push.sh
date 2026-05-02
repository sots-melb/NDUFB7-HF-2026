#!/bin/bash
PROJECT_DIR="$HOME/Projects/NDUFB7_HF_2026_04_20"
cd "$PROJECT_DIR" || exit 1
echo "========================================"
echo "V143: GitHub推送 (SSH→HTTPS回退)"
echo "========================================"
CURRENT_URL=$(git remote get-url origin 2>/dev/null || echo "")
echo "Current remote: $CURRENT_URL"
if echo "$CURRENT_URL" | grep -q "git@github.com"; then
    HTTPS_URL=$(echo "$CURRENT_URL" | sed 's|git@github.com:|https://github.com/|' | sed 's|\.git$||').git
    git remote set-url origin "$HTTPS_URL"
    echo "[DONE] Switched to HTTPS: $HTTPS_URL"
elif [ -z "$CURRENT_URL" ]; then
    git remote add origin https://github.com/sots-melb/NDUFB7-HF-2026.git 2>/dev/null || git remote set-url origin https://github.com/sots-melb/NDUFB7-HF-2026.git
    echo "[DONE] Set HTTPS remote"
else
    echo "[PASS] Remote OK"
fi
git config --local credential.helper store 2>/dev/null || true
git add -A
git commit -m "feat(20260501_v3): BIC resolution + bulk ferroptosis + Moran alternative + pub figs
- V140: 6-method BIC evidence + Methods/Discussion text
- V141: Bulk ferroptosis rho=-0.229 p=4.6e-05
- V142: Pseudo-spatial Moran I + temporal/etiological proxies
- V144: Publication-grade Fig 3 (pending)" || echo "[INFO] Nothing new to commit"
echo ""
echo "=== 推送（提示输入密码时粘贴GitHub Token）==="
git push origin main || git push origin master
echo ""
echo "[DONE] V143 complete"
