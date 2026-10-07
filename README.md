# Simulación de enlace satelital: TV, internet y voz

## Ejecutar

Requiere Python 3.10 o posterior. En una terminal, entrar a esta carpeta:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 simulador.py
```

En Windows: activar con `.venv\Scripts\activate` y utilizar `python` si corresponde.
Se abre una ventana de Matplotlib con cinco controles: atenuación de lluvia ascendente y descendente, potencia total del satélite, diámetro receptor y potencia terrestre. El botón Exportar CSV guarda el escenario vigente y un barrido de lluvia en `resultados`.

Si no aparece una ventana, ejecutar desde una terminal local con un backend gráfico de Matplotlib. En servidores o notebooks sin escritorio, usar el modo de exportación:

```bash
python3 simulador.py --exportar --salida resultados
python3 simulador.py --exportar --lluvia-u 6 --lluvia-d 5 --potencia-satelite 40 --diametro 2.4 --salida escenario_lluvia
```

Se generan `panel.png`, `resumen.csv` y `barrido_lluvia.csv`. La carpeta `ejemplo_resultados` contiene una ejecución inicial.

## Servicios simultáneos

| Servicio | Tasa bruta | Modulación | Fracción de potencia | Paquete | PER máxima | Latencia unidireccional máxima |
|---|---:|---|---:|---:|---:|---:|
| TV | 12 Mbit/s | QPSK | 55 % | 10528 bits | 1 % | 1500 ms |
| Internet | 20 Mbit/s | 16-QAM | 40 % | 12000 bits | 2 % | 600 ms |
| Voz | 0.128 Mbit/s | BPSK | 5 % | 1280 bits | 1 % | 300 ms |

Todos son valores educativos editables en `SERVICIOS`, no requisitos de estándares ni umbrales universales de calidad. Voz representa un flujo agregado con cabeceras idealizadas; TV no genera video e internet no implementa TCP. Modulación y servicio no están asociados obligatoriamente: probar otras combinaciones es parte de la actividad.

Los tres servicios ocupan portadoras separadas (FDM), comparten una potencia RF total de 100 W en tierra y 20 W en el satélite, y una banda total disponible de 36 MHz. Los pesos se normalizan; sus potencias suman exactamente el total disponible. La potencia satelital representa RF útil asignable tras cualquier reducción de operación del amplificador; no se modelan saturación ni intermodulación.

## Modelo y fórmulas

Satélite GEO transparente; dos trayectos de 38000 km, subida a 14 GHz y bajada a 12 GHz. Ganancias satelitales de 32 dBi y ganancia terrestre transmisora de 45 dBi. Ganancia receptora calculada para una parábola de 1.8 m y eficiencia 0.65.

1. `Lfs = 92.45 + 20 log10(f_GHz) + 20 log10(R_km)`.
2. `G_rx = 10 log10(eficiencia × (pi D f / c)^2)`.
3. `PIRE = P_tx_dBW - pérdida_tx + G_tx`.
4. `C_rx = PIRE - Lfs - otras_pérdidas - lluvia + G_rx`.
5. `C/N0 = C_rx - 10 log10(k T)` en dBHz. Temperaturas equivalentes: 500 K en satélite y 200 K en tierra, referidas a terminales de antena; incluyen idealmente la cadena RX. No se resta un alimentador RX adicional.
6. Para el transpondedor transparente: `1/(C/N0_total) = 1/(C/N0_u) + 1/(C/N0_d)`, en escala lineal y con las mismas referencias por portadora.
7. `Eb/N0 = C/N0_total - 10 log10(Rb)`.
8. BPSK y QPSK coherentes con Gray: `BER = 0.5 erfc(sqrt(Eb/N0_lineal))`. Para 16-QAM Gray: `BER ≈ 0.375 erfc(sqrt(0.4 Eb/N0_lineal))`; aproximación principalmente útil a BER baja.
9. `PER = 1 - (1 - BER)^L`, con L bits por paquete y errores independientes.
10. `Goodput = Rb × (1 - PER)`: tasa esperada de bits en paquetes íntegros; no descuenta cabeceras ni aplica retransmisión.
11. `B = Rb / log2(M) × (1 + rolloff)`, con rolloff 0.25. Es banda ocupada de modulación, no necesariamente ancho de ruido del receptor. Se utiliza Eb/N0 directamente para BER.
12. Margen = Eb/N0 obtenido menos Eb/N0 requerido para cumplir PER, determinado por búsqueda numérica.
13. Latencia unidireccional = `(Ru + Rd)/c + 20 ms` de procesamiento. RTT aproximado es el doble, suponiendo retorno simétrico. No incluye colas, jitter, serialización ni buffers de aplicación.

**Cumple** exige PER por debajo del objetivo, retardo por debajo del límite y banda suficiente para transmitir todas las portadoras. Si falta banda, el modelo declara el plan conjunto inviable y fija goodput en cero: no ejecuta un planificador ni redistribuye tasas.

El panel muestra tasas brutas y goodput, margen de calidad, pérdida de paquetes y sensibilidad a lluvia. El gráfico de lluvia fija para el barrido `lluvia_u = 1.2 × lluvia_d`; es una relación de escenarios asumida, no una ley física de propagación. Los controles de lluvia son independientes para el escenario actual.

## Interpretación del ejemplo

- Banda total: 13.91 MHz; latencia unidireccional: aproximadamente 273.5 ms; RTT aproximado: 547 ms.
- TV: Eb/N0 ≈ 11.01 dB; PER ≈ 0.264 %; cumple con poco margen.
- Internet: Eb/N0 ≈ 7.41 dB; PER cercana al 100 %; no cumple. Sin FEC, una BER moderada prácticamente elimina los paquetes largos íntegros.
- Voz: Eb/N0 ≈ 20.32 dB; PER muy pequeña; cumple el límite didáctico de 300 ms, aunque el retardo GEO puede afectar una conversación.

No se estima disponibilidad anual: para ello harían falta estadísticas meteorológicas locales y un modelo de propagación/tiempo. El estado Cumple representa únicamente el escenario evaluado. Tampoco se modelan códecs, FEC, ACM, TCP, jitter, interferencia, tráfico variable o fidelidad audiovisual. Los resultados son analíticos, deterministas y útiles para comparar decisiones del presupuesto de enlace.

## Actividad sugerida

1. Registrar Eb/N0, PER, margen, goodput y latencia en cielo despejado.
2. Añadir 3, 6 y 10 dB de lluvia descendente y 1.2 veces esos valores en subida.
3. Comparar duplicar potencia satelital con duplicar diámetro receptor. ¿Por qué la mejora total puede ser menor de 3 o 6 dB cuando limita la subida?
4. En `SERVICIOS`, reducir internet a 5 Mbit/s manteniendo 16-QAM; después cambiar a QPSK y comparar margen y banda.
5. Cambiar los pesos a TV 0.40, internet 0.55, voz 0.05 y explicar quién gana y quién pierde.
6. Reducir `ancho_total_mhz` a 10 para identificar un plan de banda inviable.
7. Entregar una tabla de escenarios y una propuesta justificada para los tres servicios.

## Referencias técnicas

- MathWorks, expresiones analíticas de BER para canales AWGN: https://www.mathworks.com/help/comm/ug/analytical-expressions-used-in-berawgn-function-and-bit-error-rate-analysis-app.html
- Matplotlib, controles Slider: https://matplotlib.org/stable/gallery/widgets/slider_demo.html
- Maral, Bousquet y Sun, *Satellite Communications Systems*, 6.ª edición, Wiley.
