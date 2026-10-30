# Ménsula Alféizar · repisa de ventana para gatos sin taladro ni ventosas

Diseño paramétrico y documentación del prototipo 1 (rev. A).

- `mensula.py` — generador: STL del cuerpo y de las piezas impresas, plantillas SVG 1:1 y cálculo estático. `python3 mensula.py` (requiere `trimesh shapely manifold3d numpy`).
- `salida/` — resultados del generador con los parámetros actuales:
  - `mensula_cuerpo.stl` cuerpo (contrachapado 18 mm o PETG en dos mitades)
  - `pie_pared.stl`, `mordaza_fleje.stl` piezas impresas en 3D
  - `mensula_perfil_1a1.svg` plantilla de corte del cuerpo, `fleje_gancho_1a1.svg` desarrollo del fleje de acero
  - `mensula_datos.json` parámetros y cálculo
- `documento.html` — documento de diseño: principio de funcionamiento, calculadora, rangos de ajuste, piezas con visor 3D, lista de materiales, montaje, protocolo de pruebas y medidas que faltan. Se abre en local (necesita conexión para Three.js).
- `plantilla_documento.html` — plantilla del documento sin los STL incrustados (la rellena el paso de construcción).

Principio: el gato carga un brazo en voladizo; la pata apoya en la pared y un fleje de acero con diente engancha el labio interior del marco de aluminio, que es lo que impide que deslice hacia la habitación. `N = W·L/h`.
