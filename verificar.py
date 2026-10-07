"""Verificaciones físicas y de consistencia del modelo."""
from dataclasses import replace
from simulador import Enlace, SERVICIOS, evaluar, fspl, ber_awgn
import math

e=Enlace(); base=evaluar(e)
assert math.isclose(fspl(14,76000)-fspl(14,38000),20*math.log10(2))
assert math.isclose(sum(s['potencia_fraccion'] for s in base['servicios']),1)
assert math.isclose(base['banda_usada_mhz'],13.91)
assert math.isclose(ber_awgn('BPSK',0),0.07864960352514257)
lluvia=evaluar(replace(e,lluvia_u_db=8,lluvia_d_db=8))
mejora=evaluar(replace(e,potencia_satelite_w=40,diametro_rx_m=3.6))
for b,l,m in zip(base['servicios'],lluvia['servicios'],mejora['servicios']):
    assert math.isclose(b['cn0_total_dbhz']-l['cn0_total_dbhz'],8)
    assert l['goodput_mbps'] <= b['goodput_mbps'] <= m['goodput_mbps']
    assert 0 <= b['ber'] <= .5 and 0 <= b['per'] <= 1
    assert b['cn0_total_dbhz'] < min(b['cn0_u_dbhz'],b['cn0_d_dbhz'])
    assert math.isclose(m['recibida_d_dbw']-b['recibida_d_dbw'],30*math.log10(2))
assert all(s['goodput_mbps']==0 and s['estado']=='No cumple' for s in evaluar(replace(e,ancho_total_mhz=10))['servicios'])
print('Verificaciones correctas: unidades, reparto, BER conocida, combinación de ruido, lluvia, mejoras y banda.')
