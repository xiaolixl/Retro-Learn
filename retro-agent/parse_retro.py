import json
import sys

def main():
    json_path = r"C:\Users\13558\Desktop\retro_result.json"
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    routes = data['data']['retrosynthesis_routes']
    # 按 route_rank 排序
    routes_sorted = sorted(routes, key=lambda x: x['route_rank'])
    top3 = routes_sorted[:3]
    
    descriptions = []
    for route in top3:
        rank = route['route_rank']
        score = route['score']
        # 提取反应条件
        conditions = route.get('reaction_condition', [])
        cond_text = '，'.join(conditions) if conditions else '标准条件'
        # 提取反应物 SMILES（可以简化，不展示 SMILES）
        reactants = [r['smiles'] for r in route['reactants']]
        # 可以根据需要将 SMILES 转换为名称，但这里我们只描述
        desc = f"第{rank}条路径（评分：{score}）："
        desc += f" 反应条件为 {cond_text}，"
        desc += f" 使用反应物 {', '.join(reactants)}，"
        desc += f" 通过反应模板 {route['reaction_template']} 生成目标分子。"
        descriptions.append(desc)
    
    output = "\n".join(descriptions)
    print(output)

if __name__ == '__main__':
    main()