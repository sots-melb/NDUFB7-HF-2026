import urllib.request, json
import os

url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=myocardial+infarction+time+course+human+admission&retmode=json&retmax=10"
out_file = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/STEMI_Acute_Candidates.txt")

try:
    response = urllib.request.urlopen(url)
    data = json.loads(response.read())
    ids = ",".join(data["esearchresult"]["idlist"])
    
    with open(out_file, "w", encoding="utf-8") as f:
        if ids:
            sum_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={ids}&retmode=json"
            res2 = urllib.request.urlopen(sum_url)
            sum_data = json.loads(res2.read())
            f.write("🔍 挖掘到的 STEMI 急性期候选队列:\n")
            for uid in data["esearchresult"]["idlist"]:
                info = sum_data["result"][uid]
                # 修复核心：将变量提前提取，避免在 f-string 中使用反斜杠
                acc = info.get("accession", "N/A")
                title = info.get("title", "N/A")
                f.write(f" - {acc}: {title}\n")
                print(f"✅ 找到候选队列: {acc}")
            print(f"🎉 爬取完成！结果已保存至: {out_file}")
        else:
            print("⚠️ 未找到匹配的急性期队列。")
except Exception as e:
    print(f"❌ API请求失败: {e}")
