#!/bin/bash
# V138: GitHub SSH配置 + 推送

PROJECT_DIR="$HOME/Projects/NDUFB7_HF_2026_04_20"
cd "$PROJECT_DIR" || exit 1

echo "========================================"
echo "V138: GitHub SSH配置"
echo "========================================"

# 1. 生成SSH密钥（如不存在）
if [ ! -f ~/.ssh/id_ed25519.pub ]; then
    echo "[ACTION] 生成ed25519密钥..."
    ssh-keygen -t ed25519 -C "ndufb7_hf_2026@lab" -f ~/.ssh/id_ed25519 -N ""
    echo "[DONE] 密钥已生成"
else
    echo "[PASS] SSH密钥已存在"
fi

# 2. 显示公钥（需手动复制到GitHub）
echo ""
echo "=== 公钥内容（复制到GitHub）==="
cat ~/.ssh/id_ed25519.pub
echo ""
echo "手动步骤:"
echo "1. 复制上面以 ssh-ed25519 开头的整行"
echo "2. 打开 https://github.com/settings/keys"
echo "3. 点击 'New SSH key' → 粘贴 → Title: 'NDUFB7_Lab_PC' → Add"

# 3. 启动SSH agent
eval "$(ssh-agent -s)" 2>/dev/null
ssh-add ~/.ssh/id_ed25519 2>/dev/null

# 4. 测试连接
echo ""
echo "=== 测试GitHub SSH连接 ==="
ssh -T git@github.com 2>&1 || true

# 5. 配置项目remote
if [ -d .git ]; then
    CURRENT_URL=$(git remote get-url origin 2>/dev/null || echo "")
    if echo "$CURRENT_URL" | grep -q "https://github.com/"; then
        REPO_PATH=$(echo "$CURRENT_URL" | sed 's|https://github.com/||' | sed 's|\.git$||')
        SSH_URL="git@github.com:${REPO_PATH}.git"
        git remote set-url origin "$SSH_URL"
        echo "[DONE] Remote改为SSH: $SSH_URL"
    fi
    
    # 6. 添加今晚批次并推送
    echo ""
    echo "=== 准备推送 ==="
    git add -A
    git commit -m "feat(20260501): stepwise depletion + tri-cohort modality + bulk ferroptosis + CellChat CM-CM
                    
- V120: Fig3C Monocle3 stepwise breakpoint (p=4.3e-04)
- V134: Three-cohort distribution comparison (GSE183852/121893/168742)
- V136: GSE57338 bulk ferroptosis scoring (gene-level unlock)
- V113A: CellChat CM-CM interaction (ANXA1-FPR1 dominant)
- V129/130: Results P3-4 draft + Methods supplement
- V131: Honest threshold narrative pivot (bimodal→stepwise)" || echo "[INFO] Nothing to commit or commit failed"
    
    echo ""
    echo "=== 推送命令（SSH配置完成后执行）==="
    echo "git push origin main  # 或 master"
else
    echo "[WARN] 项目目录不是Git仓库"
    echo "[ACTION] git init && git remote add origin git@github.com:YOUR_USER/NDUFB7-HF-2026.git"
fi

echo ""
echo "[DONE] V138 SSH配置完成"
echo "[NEXT] 复制公钥到GitHub后，执行: git push origin main"
