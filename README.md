# MAPF-INF404 — Implementación MaxSAT / PySAT

Repositorio con la implementación en Python (PySAT) de la codificación MaxSAT para Multi-Agent Path Finding (MAPF) inspirada en Asín et al. (JAIR 2022).

Trabajo realizado por Matías Sandoval y Fernando Salgado.

Resumen
-------
Este repositorio contiene código y notebooks para generar las fórmulas CNF/WCNF que modelan MAPF bajo la codificación estudiada, buscar el makespan mínimo mediante SAT (Fase 1) y optimizar la suma de costes con MaxSAT (Fase 2, usando RC2 de PySAT). También incluye utilidades para visualizar de cada isntancia.

Dependencias
------------
Instalar las librerías depende de usted y del gestor de paquetes que prefiera (pip, conda, etc.). Las piezas mínimas necesarias son:

- `pysat`

Cómo ejecutar
-------------
1. Abra `mapf.ipynb` con Jupyter Notebook o JupyterLab.
2. Ejecute las celdas en el orden en que aparecen (Run → Run All). El notebook asume que las definiciones previas ya fueron ejecutadas cuando se lanzan los experimentos.
