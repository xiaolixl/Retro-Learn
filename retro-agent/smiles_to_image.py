#!/usr/bin/env python3
"""
使用RDKit将SMILES转换为化合物结构图片。
支持单个SMILES或JSON文件（包含SMILES列表）。
"""

import sys
import os
import json
import argparse
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

def smiles_to_image(smiles, output_path, img_size=(300, 200)):
    """
    将单个SMILES字符串转换为图片。
    返回 (success, message)
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, f"无效的SMILES: {smiles}"
        # 绘制分子结构
        img = Draw.MolToImage(mol, size=img_size)
        # 确保输出目录存在
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        img.save(output_path)
        return True, output_path
    except Exception as e:
        return False, f"转换失败: {e}"

def process_smiles_list(smiles_list, output_dir, prefix="mol"):
    """
    处理SMILES列表，为每个生成图片。
    返回列表，每个元素为 (smiles, success, image_path_or_error)
    """
    results = []
    for i, smi in enumerate(smiles_list, 1):
        # 生成文件名
        # 使用简单编号，也可以基于SMILES哈希
        filename = f"{prefix}_{i}.png"
        output_path = os.path.join(output_dir, filename)
        success, msg = smiles_to_image(smi, output_path)
        results.append({
            "smiles": smi,
            "success": success,
            "path": output_path if success else None,
            "error": None if success else msg
        })
    return results

def main():
    parser = argparse.ArgumentParser(description="将SMILES转换为化合物结构图片")
    parser.add_argument("-s", "--smiles", help="单个SMILES字符串")
    parser.add_argument("-j", "--json", help="包含SMILES列表的JSON文件路径")
    parser.add_argument("-o", "--output", default="./smiles_images", help="输出目录（默认: ./smiles_images）")
    parser.add_argument("--prefix", default="mol", help="输出图片文件名前缀（默认: mol）")
    args = parser.parse_args()
    
    if not args.smiles and not args.json:
        parser.print_help()
        sys.exit(1)
    
    # 确保输出目录存在
    os.makedirs(args.output, exist_ok=True)
    
    if args.smiles:
        # 单个SMILES
        output_path = os.path.join(args.output, f"{args.prefix}_1.png")
        success, msg = smiles_to_image(args.smiles, output_path)
        if success:
            print(json.dumps({"success": True, "path": output_path}))
        else:
            print(json.dumps({"success": False, "error": msg}))
    elif args.json:
        # 从JSON文件读取SMILES列表
        if not os.path.exists(args.json):
            print(json.dumps({"success": False, "error": f"JSON文件不存在: {args.json}"}))
            sys.exit(1)
        with open(args.json, 'r') as f:
            data = json.load(f)
        # 假设JSON是一个列表，或者有'reactants'字段
        smiles_list = []
        if isinstance(data, list):
            # 如果列表元素是字符串，直接使用；如果是字典，可能需要提取smiles字段
            for item in data:
                if isinstance(item, str):
                    smiles_list.append(item)
                elif isinstance(item, dict) and 'smiles' in item:
                    smiles_list.append(item['smiles'])
                else:
                    # 尝试其他常见字段
                    pass
        elif isinstance(data, dict):
            # 尝试常见字段
            if 'reactants' in data:
                reactants = data['reactants']
                if isinstance(reactants, list):
                    smiles_list.extend(reactants)
            elif 'smiles' in data:
                smiles_list.append(data['smiles'])
        if not smiles_list:
            print(json.dumps({"success": False, "error": "未找到有效的SMILES列表"}))
            sys.exit(1)
        # 处理列表
        results = process_smiles_list(smiles_list, args.output, args.prefix)
        print(json.dumps({"success": True, "results": results}))

if __name__ == "__main__":
    main()