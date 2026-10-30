#!/usr/bin/env python3
"""
Ménsula de alféizar para gatos · generador paramétrico
======================================================
Genera, a partir de unos pocos parámetros en milímetros:
  - STL del cuerpo de la ménsula (contrachapado 18 mm cortado con plantilla o CNC; también imprimible en dos mitades)
  - STL de las piezas impresas en 3D: pie de pared regulable y mordaza del fleje
  - SVG 1:1 del perfil del cuerpo (plantilla de corte) y del desarrollo del fleje-gancho de acero
  - JSON con los parámetros y el cálculo estático (fuerzas, tensiones, márgenes)

Sistema de coordenadas del perfil (vista lateral):
  x → hacia la habitación, y → arriba, origen en la arista interior del alféizar (donde el alféizar toca la pared).
  El cuerpo se extruye en z (espesor de la tabla).

Uso:  python3 mensula.py [carpeta_salida]
Requiere: trimesh, shapely, manifold3d, numpy
"""
import json, math, sys
from pathlib import Path
import numpy as np
import trimesh
from shapely.geometry import Polygon, box
from shapely.ops import unary_union

# ---------------------------------------------------------------------------
# 1. PARÁMETROS (mm, kgf). Cambia aquí y vuelve a ejecutar.
# ---------------------------------------------------------------------------
P = dict(
    # Cuerpo (contrachapado de abedul 18 mm o PETG impreso)
    espesor=18.0,            # espesor del cuerpo (z)
    brazo_largo=320.0,       # voladizo del brazo desde la pared hasta la punta
    brazo_alto=43.0,         # altura total del brazo (la tabla apoya encima)
    rebaje_fleje=3.0,        # rebaje superior donde va el fleje de acero (x <= rebaje_hasta)
    rebaje_hasta=100.0,
    lengueta=30.0,           # cuánto entra el cuerpo sobre el alféizar (el fleje añade el resto)
    holgura_alfeizar=3.0,    # goma bajo la lengüeta (apoyo sobre el alféizar)
    pata_alto=320.0,         # altura de la pata que baja por la pared (h del par de fuerzas)
    pata_ancho=40.0,
    pata_separacion=10.0,    # separación pata-pared (solo toca el pie regulable)
    tirante_ancho=35.0,      # diagonal de rigidez
    radio_redondeo=4.0,
    # Fleje-gancho de acero (pletina 50x3 mm, acero al carbono, doblada a 90°)
    fleje_ancho=50.0, fleje_espesor=3.0, diente=25.0, fleje_ranura=100.0, fleje_recorrido=190.0,
    tornillos_fleje_x=(20.0, 80.0), tornillo_fleje_d=6.5,     # M6 pasantes en vertical
    # Pie de pared regulable (M8 + taco de goma)
    pie_diametro=44.0, pie_alto=14.0, pie_y=-300.0, tornillo_pie_d=8.5,
    # Tabla
    tabla_fondo=320.0, tabla_largo=600.0, tabla_espesor=15.0, separacion_mensulas=450.0,
    # Carga de diseño
    peso_gato=6.0, factor_dinamico=3.0, coef_seguridad=1.6,
    # Perfil de aluminio de la ventana (para comprobar el labio)
    labio_espesor=1.4, labio_brazo=8.0,
)
P['carga_diseno'] = P['peso_gato'] * P['factor_dinamico'] * P['coef_seguridad']   # ≈ 29 kgf

G = 9.81

# ---------------------------------------------------------------------------
# 2. CÁLCULO ESTÁTICO
# ---------------------------------------------------------------------------
def calculo(p=P, L=None, h=None, W=None):
    """Equilibrio del par de ménsulas. L: posición de la carga desde la pared; h: altura del par (pie); W: carga total (kgf)."""
    L = p['brazo_largo'] if L is None else L
    h = -p['pie_y'] if h is None else h
    W = p['carga_diseno'] if W is None else W
    N = W * L / h                     # empuje total de los pies contra la pared = tiro horizontal del gancho (kgf)
    por_mensula = dict(W=W / 2, N=N / 2)
    F = por_mensula['N'] * G          # N por gancho
    # Diente del fleje: ménsula corta de longitud ~diente/2 con carga F, sección b·t
    Z_fleje = p['fleje_ancho'] * p['fleje_espesor'] ** 2 / 6
    sigma_fleje = F * (p['diente'] * 0.6) / Z_fleje                 # MPa
    # Labio de aluminio del marco: voladizo labio_brazo con carga F repartida en el ancho del fleje
    Z_labio = p['fleje_ancho'] * p['labio_espesor'] ** 2 / 6
    sigma_labio = F * p['labio_brazo'] / Z_labio                    # MPa (6063-T5 cede ≈ 145–190 MPa)
    # Tornillos M6 del fleje a cortante (2 por ménsula)
    tau_tornillo = F / (2 * math.pi * 3.0 ** 2)                      # MPa (acero 8.8 admite > 200)
    # Presión del pie sobre la pared
    area_pie = math.pi * (p['pie_diametro'] / 2) ** 2 / 100          # cm²
    presion_pie = por_mensula['N'] / area_pie                        # kgf/cm² (yeso/drywall tolera > 5)
    # Brazo de contrachapado como voladizo simple (sin contar el tirante, lado seguro)
    M_brazo = por_mensula['W'] * G * L                                # N·mm
    Z_brazo = p['espesor'] * p['brazo_alto'] ** 2 / 6
    sigma_brazo = M_brazo / Z_brazo                                   # MPa (abedul ≈ 40 MPa admisible)
    # Fricción útil bajo la lengüeta (goma sobre piedra/aluminio, μ≈0.6) que descarga al gancho
    friccion = 0.6 * W
    return dict(L=L, h=h, W=W, N=N, por_mensula=por_mensula, F_gancho_N=F,
                sigma_fleje_MPa=sigma_fleje, sigma_labio_MPa=sigma_labio, tau_tornillo_MPa=tau_tornillo,
                presion_pie_kgf_cm2=presion_pie, sigma_brazo_MPa=sigma_brazo, friccion_lengueta_kgf=friccion,
                N_neto_kgf=max(0.0, N - friccion),
                ok=dict(fleje=sigma_fleje < 150, labio=sigma_labio < 100, tornillo=tau_tornillo < 100,
                        pie=presion_pie < 3, brazo=sigma_brazo < 30))

# ---------------------------------------------------------------------------
# 3. PERFIL 2D DEL CUERPO
# ---------------------------------------------------------------------------
def perfil_cuerpo(p=P):
    y0 = p['holgura_alfeizar']
    brazo = box(-p['lengueta'], y0, p['brazo_largo'], p['brazo_alto'])
    rebaje = box(-p['lengueta'] - 1, p['brazo_alto'] - p['rebaje_fleje'], p['rebaje_hasta'], p['brazo_alto'] + 1)
    pata = box(p['pata_separacion'], -p['pata_alto'], p['pata_separacion'] + p['pata_ancho'], y0 + 1)
    # Tirante diagonal desde el brazo (x≈0.9·L) hasta la pata (y≈-0.85·h)
    a = np.array([p['brazo_largo'] * 0.9, y0 + 1]); b = np.array([p['pata_separacion'] + p['pata_ancho'] - 1, -p['pata_alto'] * 0.85])
    d = (b - a); d /= np.linalg.norm(d); n = np.array([-d[1], d[0]]) * p['tirante_ancho'] / 2
    tirante = Polygon([a + n, b + n, b - n, a - n])
    cuerpo = unary_union([brazo, pata, tirante]).difference(rebaje)
    r = p['radio_redondeo']
    cuerpo = cuerpo.buffer(r, join_style=1).buffer(-2 * r, join_style=1).buffer(r, join_style=1)  # redondea cóncavos y convexos
    return cuerpo

def cilindro_eje(eje, d, largo, centro):
    """Cilindro de diámetro d a lo largo del eje 'x' o 'y', centrado en 'centro' (x, y, z)."""
    c = trimesh.creation.cylinder(radius=d / 2, height=largo, sections=48)
    if eje == 'y': c.apply_transform(trimesh.transformations.rotation_matrix(math.pi / 2, [1, 0, 0]))
    if eje == 'x': c.apply_transform(trimesh.transformations.rotation_matrix(math.pi / 2, [0, 1, 0]))
    c.apply_translation(centro); return c

def cuerpo_3d(p=P):
    cuerpo = trimesh.creation.extrude_polygon(perfil_cuerpo(p), p['espesor'])
    zc = p['espesor'] / 2
    taladros = [cilindro_eje('y', p['tornillo_fleje_d'], p['brazo_alto'] * 2, (x, p['brazo_alto'] / 2, zc)) for x in p['tornillos_fleje_x']]
    taladros.append(cilindro_eje('x', p['tornillo_pie_d'], p['pata_ancho'] * 3, (p['pata_separacion'] + p['pata_ancho'] / 2, p['pie_y'], zc)))
    return trimesh.boolean.difference([cuerpo] + taladros, engine='manifold')

def pie_pared(p=P):
    """Pie regulable: cilindro con alojamiento hexagonal para cabeza M8 por detrás y cazoleta para disco de goma por delante."""
    cuerpo = trimesh.creation.cylinder(radius=p['pie_diametro'] / 2, height=p['pie_alto'], sections=96)
    hexagono = trimesh.creation.cylinder(radius=13.2 / math.sqrt(3), height=7, sections=6)   # 13 mm entre caras + holgura
    hexagono.apply_translation((0, 0, -p['pie_alto'] / 2 + 3.5 - 0.01))
    pasante = trimesh.creation.cylinder(radius=4.4, height=p['pie_alto'] * 2, sections=48)
    cazoleta = trimesh.creation.cylinder(radius=p['pie_diametro'] / 2 - 2, height=3, sections=96)
    cazoleta.apply_translation((0, 0, p['pie_alto'] / 2 - 1.5 + 0.01))
    return trimesh.boolean.difference([cuerpo, hexagono, pasante, cazoleta], engine='manifold')

def mordaza(p=P):
    """Barra de presión que reparte el apriete de las palomillas M6 sobre la ranura del fleje."""
    largo = p['tornillos_fleje_x'][1] - p['tornillos_fleje_x'][0] + 40
    barra = trimesh.creation.box((largo, p['fleje_ancho'] - 4, 6))
    taladros = []
    for dx in (-(p['tornillos_fleje_x'][1] - p['tornillos_fleje_x'][0]) / 2, (p['tornillos_fleje_x'][1] - p['tornillos_fleje_x'][0]) / 2):
        c = trimesh.creation.cylinder(radius=p['tornillo_fleje_d'] / 2, height=20, sections=48); c.apply_translation((dx, 0, 0)); taladros.append(c)
    return trimesh.boolean.difference([barra] + taladros, engine='manifold')

# ---------------------------------------------------------------------------
# 4. DIBUJOS SVG (1:1 en mm; 1 unidad SVG = 1 mm)
# ---------------------------------------------------------------------------
def svg_cabecera(w, h, titulo):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}mm" height="{h}mm" viewBox="0 0 {w} {h}" font-family="sans-serif">'
            f'<title>{titulo}</title><rect width="{w}" height="{h}" fill="white"/>')

def svg_perfil(p=P):
    poly = perfil_cuerpo(p)
    minx, miny, maxx, maxy = poly.bounds
    m = 15; w = maxx - minx + 2 * m; h = maxy - miny + 2 * m
    # Volteamos y (SVG crece hacia abajo)
    tx = lambda x: x - minx + m; ty = lambda y: maxy - y + m
    pts = ' '.join(f'{tx(x):.2f},{ty(y):.2f}' for x, y in poly.exterior.coords)
    s = [svg_cabecera(w, h, 'Perfil de la ménsula · plantilla 1:1')]
    s.append(f'<polygon points="{pts}" fill="none" stroke="black" stroke-width="0.4"/>')
    for x in p['tornillos_fleje_x']:
        s.append(f'<line x1="{tx(x):.2f}" y1="{ty(p["brazo_alto"]):.2f}" x2="{tx(x):.2f}" y2="{ty(p["holgura_alfeizar"]):.2f}" stroke="black" stroke-width="0.3" stroke-dasharray="2 1.5"/>')
        s.append(f'<text x="{tx(x) + 2:.1f}" y="{ty(p["brazo_alto"]) - 2:.1f}" font-size="4">Ø{p["tornillo_fleje_d"]} vertical (M6)</text>')
    cx, cy = tx(p['pata_separacion'] + p['pata_ancho'] / 2), ty(p['pie_y'])
    s.append(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{p["tornillo_pie_d"] / 2:.2f}" fill="none" stroke="black" stroke-width="0.3"/>')
    s.append(f'<text x="{cx + 8:.1f}" y="{cy + 1.5:.1f}" font-size="4">Ø{p["tornillo_pie_d"]} horizontal (tuerca de golpe M8, pie)</text>')
    # Referencias
    s.append(f'<line x1="{tx(0):.2f}" y1="{ty(maxy) - 5:.2f}" x2="{tx(0):.2f}" y2="{ty(miny) + 5:.2f}" stroke="#888" stroke-width="0.3" stroke-dasharray="4 2"/>')
    s.append(f'<text x="{tx(0) + 2:.1f}" y="{ty(miny) + 4:.1f}" font-size="4" fill="#555">plano de la pared</text>')
    s.append(f'<line x1="{tx(-p["lengueta"]) - 5:.2f}" y1="{ty(0):.2f}" x2="{tx(0):.2f}" y2="{ty(0):.2f}" stroke="#888" stroke-width="0.3" stroke-dasharray="4 2"/>')
    s.append(f'<text x="{tx(-p["lengueta"]) - 4:.1f}" y="{ty(0) + 6:.1f}" font-size="4" fill="#555">alféizar</text>')
    s.append(f'<text x="{m}" y="{h - 4:.1f}" font-size="5">Ménsula de alféizar · perfil 1:1 · contrachapado {p["espesor"]:.0f} mm · brazo {p["brazo_largo"]:.0f} · pata {p["pata_alto"]:.0f} · imprimir al 100 % (sin ajustar a página)</text>')
    s.append('</svg>'); return '\n'.join(s)

def svg_fleje(p=P):
    """Desarrollo plano del fleje-gancho: tramo recto + diente, con la ranura de regulación. Doblar 90° por la línea."""
    L = p['fleje_recorrido'] + p['diente'] + p['fleje_espesor']  # tolerancia de doblado ≈ espesor
    w = p['fleje_ancho']; m = 15
    s = [svg_cabecera(L + 2 * m, w + 2 * m + 20, 'Fleje-gancho · desarrollo 1:1')]
    s.append(f'<rect x="{m}" y="{m}" width="{L:.1f}" height="{w}" fill="none" stroke="black" stroke-width="0.4"/>')
    xr = m + 30; s.append(f'<rect x="{xr}" y="{m + w / 2 - p["tornillo_fleje_d"] / 2:.2f}" width="{p["fleje_ranura"]}" height="{p["tornillo_fleje_d"]}" rx="{p["tornillo_fleje_d"] / 2}" fill="none" stroke="black" stroke-width="0.4"/>')
    xd = m + p['fleje_recorrido']
    s.append(f'<line x1="{xd:.1f}" y1="{m - 5}" x2="{xd:.1f}" y2="{m + w + 5}" stroke="black" stroke-width="0.3" stroke-dasharray="3 2"/>')
    s.append(f'<text x="{xd - 18:.1f}" y="{m + w + 11}" font-size="4">doblar 90° hacia abajo</text>')
    s.append(f'<text x="{xd + 2:.1f}" y="{m - 2}" font-size="4">diente {p["diente"]:.0f}</text>')
    s.append(f'<text x="{m}" y="{m + w + 18}" font-size="5">Pletina acero {p["fleje_ancho"]:.0f}×{p["fleje_espesor"]:.0f} mm · largo total {L:.0f} mm · ranura {p["fleje_ranura"]:.0f}×{p["tornillo_fleje_d"]} · 2 unidades · desbarbar y pintar o galvanizar</text>')
    s.append('</svg>'); return '\n'.join(s)

# ---------------------------------------------------------------------------
# 5. SALIDA
# ---------------------------------------------------------------------------
def main(salida='salida'):
    out = Path(salida); out.mkdir(parents=True, exist_ok=True)
    piezas = {'mensula_cuerpo.stl': cuerpo_3d(), 'pie_pared.stl': pie_pared(), 'mordaza_fleje.stl': mordaza()}
    informe = {}
    for nombre, malla in piezas.items():
        assert malla.is_watertight, f'{nombre} no es estanca'
        malla.export(out / nombre)
        informe[nombre] = dict(volumen_cm3=round(malla.volume / 1000, 1), caja_mm=[round(v, 1) for v in malla.extents])
    (out / 'mensula_perfil_1a1.svg').write_text(svg_perfil(), encoding='utf-8')
    (out / 'fleje_gancho_1a1.svg').write_text(svg_fleje(), encoding='utf-8')
    datos = dict(parametros=P, calculo=calculo(), piezas=informe,
                 calculo_tabla=[calculo(L=L, h=h) for L in (250, 320, 350) for h in (250, 320, 400)])
    (out / 'mensula_datos.json').write_text(json.dumps(datos, indent=1, ensure_ascii=False), encoding='utf-8')
    c = datos['calculo']
    print(f"Carga de diseño {c['W']:.1f} kgf · empuje pies/tiro gancho N = {c['N']:.1f} kgf ({c['por_mensula']['N']:.1f} por ménsula)")
    print(f"σ fleje {c['sigma_fleje_MPa']:.0f} MPa · σ labio aluminio {c['sigma_labio_MPa']:.0f} MPa · τ tornillo {c['tau_tornillo_MPa']:.1f} MPa · presión pie {c['presion_pie_kgf_cm2']:.2f} kgf/cm² · σ brazo {c['sigma_brazo_MPa']:.1f} MPa")
    print('Comprobaciones:', c['ok']); print('Piezas:', json.dumps(informe, ensure_ascii=False))

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).parent / 'salida'))
