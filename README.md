# 🌈 RGBC-SpectraSim

**A Python-based Interactive RGBC Display Gamut & Metamerism Simulator**

[![Python 3.x](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

RGBC-SpectraSim 是一个轻量级的桌面端交互式模拟器，用于探索四基色（Red, Green, Cyan, Blue）显示系统（如先进的量子点或微型 LED 面板）的色彩映射逻辑。

在传统的 RGB 三基色系统中，色域内的坐标与驱动电平是一一对应的。但引入第四基色（Cyan）后，系统成为了**超定系统（Underdetermined System）**，即同一个 CIE 1931 坐标可以由无数种物理光谱组合而成。本项目旨在通过可视化手段，拆解并分析不同的色域映射算法（GMA）对底层光谱分布的实际影响。

## ✨ 核心特性 (Features)

* **精准的色彩空间约束**：内置精简版 CIE 1931 光谱轨迹，利用重心坐标（Barycentric Coordinates）算法实现物理级边界碰撞检测，确保交互点始终被限制在合法的 RGBC 四边形色域内。
* **实时窄带光谱渲染**：模拟现代高色域面板特性，生成半峰宽极窄的高斯光谱（$\sigma \approx 8-10nm$）。
* **动态彩虹渐变 SPD**：基于物理波长估算算法，在光谱辐射功率（SPD）图表底层实时渲染对应的可见光色彩渐变，直观展示能量分布。
* **三大底层驱动策略 (GMA Modes)**：
  * 🗡️ **革新派 (R-C Split)**：最大化色纯度，沿红-青对角线切割，优先使用 Cyan 灯珠合成中心白点。
  * 🛡️ **保守派 (G-B Split)**：沿绿-蓝对角线切割，保留传统 RGB 调色逻辑。
  * 🌍 **环保少女How dare you派 (4-Color Blend)**：四色全开，均摊单颗 LED 功耗，实现面板全局亮度的最大化增强（Brightness Enhancement）。

## 🚀 安装与运行 (Installation & Usage)

依赖纯粹的科学计算与绘图基础库。

1. **克隆仓库**：
   ```bash
   git clone [https://github.com/yourusername/RGBC-SpectraSim.git](https://github.com/yourusername/RGBC-SpectraSim.git)
   cd RGBC-SpectraSim
