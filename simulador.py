"""Enlace GEO transparente multicarrier: modelo analítico educativo."""
from dataclasses import dataclass, replace
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np

@dataclass
class Enlace:
    distancia_km: float = 38000
    frecuencia_u_ghz: float = 14
    frecuencia_d_ghz: float = 12
    potencia_tierra_w: float = 100
    potencia_satelite_w: float = 20
    ganancia_tierra_tx_dbi: float = 45
    ganancia_satelite_rx_dbi: float = 32
    ganancia_satelite_tx_dbi: float = 32
    diametro_rx_m: float = 1.8
    eficiencia_rx: float = 0.65
    temperatura_u_k: float = 500
    temperatura_d_k: float = 200
    perdida_tx_u_db: float = 2
    perdida_tx_d_db: float = 1
    otras_u_db: float = 2
    otras_d_db: float = 2
    lluvia_u_db: float = 0
    lluvia_d_db: float = 0
    ancho_total_mhz: float = 36
    rolloff: float = 0.25
    procesamiento_ms: float = 20

@dataclass
class Servicio:
    nombre: str
    tasa_mbps: float
    modulacion: str
    peso_potencia: float
    paquete_bits: int
    per_max: float
    latencia_max_ms: float

SERVICIOS = [Servicio('TV',12,'QPSK',0.55,10528,0.01,1500),
             Servicio('Internet',20,'16QAM',0.40,12000,0.02,600),
             Servicio('Voz',0.128,'BPSK',0.05,1280,0.01,300)]
COLORES = ['#2563eb','#e58916','#119977']

def db(x):
    return 10 * math.log10(x)

def fspl(f_ghz, r_km):
    return 92.45 + 20*math.log10(f_ghz) + 20*math.log10(r_km)

def ber_awgn(modulacion, ebn0_db):
    gamma = 10**(ebn0_db/10)
    if modulacion in ('BPSK','QPSK'):
        return 0.5*math.erfc(math.sqrt(gamma))
    if modulacion == '16QAM':
        # Aproximación Gray, más precisa en la región de BER baja.
        return 0.375*math.erfc(math.sqrt(0.4*gamma))
    raise ValueError('Modulación no soportada')

def validar(e, servicios):
    positivos = (e.distancia_km,e.frecuencia_u_ghz,e.frecuencia_d_ghz,
                 e.potencia_tierra_w,e.potencia_satelite_w,e.diametro_rx_m,
                 e.temperatura_u_k,e.temperatura_d_k,e.ancho_total_mhz)
    if any(v <= 0 for v in positivos) or not 0 < e.eficiencia_rx <= 1:
        raise ValueError('Distancia, frecuencias, potencias, temperaturas y tamaños deben ser positivos; 0 < eficiencia <= 1.')
    if not servicios or any(s.tasa_mbps <= 0 or s.peso_potencia <= 0 or s.paquete_bits <= 0 for s in servicios):
        raise ValueError('Cada servicio requiere tasa, peso y tamaño de paquete positivos.')
    if e.rolloff < 0 or e.procesamiento_ms < 0 or min(e.lluvia_u_db,e.lluvia_d_db) < 0:
        raise ValueError('Rolloff, retardos y pérdidas de lluvia no pueden ser negativos.')
    if any(not 0 < s.per_max < 1 or s.latencia_max_ms <= 0 for s in servicios):
        raise ValueError('Los objetivos PER y latencia deben ser válidos.')

def requerido_ebn0(s):
    # Umbral derivado del objetivo PER, sin inventar ganancia de codificación.
    lo, hi = -30., 40.
    for _ in range(80):
        mid = (lo+hi)/2
        ber = ber_awgn(s.modulacion,mid)
        per = -math.expm1(s.paquete_bits*math.log1p(-ber))
        if per > s.per_max: lo = mid
        else: hi = mid
    return hi

def evaluar(e, servicios=SERVICIOS):
    validar(e, servicios)
    k_db = db(1.380649e-23)
    gr = db(e.eficiencia_rx*(math.pi*e.diametro_rx_m*e.frecuencia_d_ghz*1e9/299792458)**2)
    lu, ld = fspl(e.frecuencia_u_ghz,e.distancia_km), fspl(e.frecuencia_d_ghz,e.distancia_km)
    total_peso = sum(s.peso_potencia for s in servicios)
    bits = {'BPSK':1,'QPSK':2,'16QAM':4}
    bw_total = sum(s.tasa_mbps*1e6/bits[s.modulacion]*(1+e.rolloff) for s in servicios)
    bw_ok = bw_total <= e.ancho_total_mhz*1e6
    retardo = 2*e.distancia_km*1000/299792458*1000 + e.procesamiento_ms
    filas=[]
    for s in servicios:
        fraccion = s.peso_potencia/total_peso
        pu,pd = db(e.potencia_tierra_w*fraccion),db(e.potencia_satelite_w*fraccion)
        eirpu = pu-e.perdida_tx_u_db+e.ganancia_tierra_tx_dbi
        eirpd = pd-e.perdida_tx_d_db+e.ganancia_satelite_tx_dbi
        # Referencia en terminales de antena RX; temperatura equivalente en ese punto.
        cu = eirpu-lu-e.otras_u_db-e.lluvia_u_db+e.ganancia_satelite_rx_dbi
        cd = eirpd-ld-e.otras_d_db-e.lluvia_d_db+gr
        cn0u = cu-k_db-db(e.temperatura_u_k)
        cn0d = cd-k_db-db(e.temperatura_d_k)
        cn0 = -db(10**(-cn0u/10)+10**(-cn0d/10))
        rb = s.tasa_mbps*1e6
        bw = rb/bits[s.modulacion]*(1+e.rolloff)
        eb = cn0-db(rb)
        ber = ber_awgn(s.modulacion,eb)
        per = -math.expm1(s.paquete_bits*math.log1p(-ber))
        margen = eb-requerido_ebn0(s)
        estado = 'Cumple' if bw_ok and per <= s.per_max and retardo <= s.latencia_max_ms else 'No cumple'
        filas.append(dict(servicio=s.nombre,modulacion=s.modulacion,potencia_fraccion=fraccion,
            tasa_mbps=s.tasa_mbps,banda_mhz=bw/1e6,pire_u_dbw=eirpu,pire_d_dbw=eirpd,
            perdida_u_db=lu,perdida_d_db=ld,recibida_u_dbw=cu,recibida_d_dbw=cd,
            cn0_u_dbhz=cn0u,cn0_d_dbhz=cn0d,cn0_total_dbhz=cn0,ebn0_db=eb,
            ber=ber,per=per,margen_db=margen,goodput_mbps=s.tasa_mbps*(1-per) if bw_ok else 0,
            latencia_ms=retardo,per_objetivo=s.per_max,latencia_objetivo_ms=s.latencia_max_ms,
            estado=estado))
    return dict(servicios=filas,banda_usada_mhz=bw_total/1e6,banda_disponible_mhz=e.ancho_total_mhz,
                banda_ok=bw_ok,ganancia_rx_dbi=gr,latencia_ms=retardo)

def exportar(e, carpeta):
    carpeta=Path(carpeta); carpeta.mkdir(parents=True,exist_ok=True)
    r=evaluar(e)
    with (carpeta/'resumen.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(r['servicios'][0])); w.writeheader(); w.writerows(r['servicios'])
    with (carpeta/'barrido_lluvia.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['lluvia_d_db','lluvia_u_db','servicio','ebn0_db','ber','per','goodput_mbps','margen_db','estado']); w.writeheader()
        for lluvia in np.linspace(0,20,101):
            for s in evaluar(replace(e,lluvia_d_db=float(lluvia),lluvia_u_db=1.2*float(lluvia)))['servicios']:
                w.writerow(dict(lluvia_d_db=lluvia,lluvia_u_db=1.2*lluvia,**{k:s[k] for k in ['servicio','ebn0_db','ber','per','goodput_mbps','margen_db','estado']}))
    return r

def figura(e, interactiva=True):
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button
    fig, axs=plt.subplots(2,2,figsize=(15,9))
    fig.subplots_adjust(left=.07,right=.97,top=.86,bottom=.32 if interactiva else .12,hspace=.50,wspace=.25)
    titulo=fig.suptitle('Enlace GEO · TV / Internet / Voz',fontsize=19,weight='bold')
    resumen=fig.text(.07,.90,'',fontsize=10)
    pie=fig.text(.07,.275 if interactiva else .02,'',fontsize=9,family='monospace')
    def dibujar(actual):
        r=evaluar(actual); ss=r['servicios']; nombres=[s['servicio'] for s in ss]
        for ax in axs.flat: ax.clear(); ax.grid(alpha=.2,axis='y')
        ax=axs[0,0]; x=np.arange(3)
        ax.bar(x-.18,[s['tasa_mbps'] for s in ss],.36,label='Tasa bruta',color='#cbd5e1')
        ax.bar(x+.18,[s['goodput_mbps'] for s in ss],.36,label='Paquetes íntegros',color=COLORES)
        ax.set_xticks(x,nombres); ax.set_ylabel('Mbit/s'); ax.set_title('Goodput esperado sin retransmisión'); ax.legend(fontsize=8)
        ax=axs[0,1]; ax.bar(nombres,[s['margen_db'] for s in ss],color=COLORES); ax.axhline(0,color='red',ls='--'); ax.set_ylabel('dB'); ax.set_title('Margen sobre el objetivo de pérdida de paquetes')
        ax=axs[1,0]
        for i,s in enumerate(ss):
            ax.bar(i,max(s['per'],1e-12),color=COLORES[i]); ax.plot([i-.3,i+.3],[s['per_objetivo']]*2,'k--')
        ax.set_yscale('log'); ax.set_ylim(1e-12,1.2); ax.set_xticks(x,nombres); ax.set_ylabel('Probabilidad PER'); ax.set_title('Pérdida de paquetes (piso visual 10⁻¹²)')
        ax=axs[1,1]; lluvias=np.linspace(0,20,81); curvas=[[] for _ in ss]
        for l in lluvias:
            rr=evaluar(replace(actual,lluvia_d_db=float(l),lluvia_u_db=1.2*float(l)))
            for i,s in enumerate(rr['servicios']): curvas[i].append(100*s['goodput_mbps']/s['tasa_mbps'])
        for i in range(3): ax.plot(lluvias,curvas[i],color=COLORES[i],label=nombres[i])
        ax.set_ylim(-2,102); ax.set_xlabel('Atenuación por lluvia descendente (dB)'); ax.set_ylabel('Goodput / tasa bruta (%)'); ax.set_title('Barrido: lluvia ascendente = 1.2 × descendente'); ax.legend(fontsize=8)
        resumen.set_text(f"Banda usada: {r['banda_usada_mhz']:.2f} / {r['banda_disponible_mhz']:.1f} MHz | RX: {r['ganancia_rx_dbi']:.1f} dBi | Latencia unidireccional: {r['latencia_ms']:.1f} ms | RTT estimado: {2*r['latencia_ms']:.1f} ms")
        pie.set_text('\n'.join(f"{s['servicio']:8} {s['modulacion']:5} | Eb/N0 {s['ebn0_db']:5.1f} dB | BER {s['ber']:.2e} | PER {s['per']:.2e} | {s['estado']}" for s in ss))
        titulo.set_text('Enlace GEO · TV / Internet / Voz'+(' · BANDA INSUFICIENTE' if not r['banda_ok'] else ''))
        fig.canvas.draw_idle()
    dibujar(e)
    if interactiva:
        especificaciones=[('Lluvia subida (dB)',0,25,e.lluvia_u_db,'lluvia_u_db'),('Lluvia bajada (dB)',0,20,e.lluvia_d_db,'lluvia_d_db'),('Potencia satélite (W)',1,80,e.potencia_satelite_w,'potencia_satelite_w'),('Diámetro receptor (m)',.4,4,e.diametro_rx_m,'diametro_rx_m'),('Potencia tierra (W)',10,500,e.potencia_tierra_w,'potencia_tierra_w')]
        sliders=[]
        for i,(label,lo,hi,val,campo) in enumerate(especificaciones):
            sliders.append((campo,Slider(fig.add_axes([.24,.205-i*.032,.58,.018]),label,lo,hi,valinit=val)))
        def actual(): return replace(e,**{campo:float(sl.val) for campo,sl in sliders})
        for _,sl in sliders: sl.on_changed(lambda _: dibujar(actual()))
        boton=Button(fig.add_axes([.85,.05,.12,.055]),'Exportar CSV')
        boton.on_clicked(lambda _: (exportar(actual(),'resultados'),print('Resultados exportados en ./resultados')))
        # Mantener controles vivos durante la sesión.
        fig._controles=(sliders,boton)
    return fig

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exportar',action='store_true',help='Generar CSV y PNG sin ventana')
    p.add_argument('--salida',default='resultados')
    p.add_argument('--lluvia-u',type=float,default=0)
    p.add_argument('--lluvia-d',type=float,default=0)
    p.add_argument('--potencia-satelite',type=float,default=20)
    p.add_argument('--diametro',type=float,default=1.8)
    a=p.parse_args()
    e=Enlace(lluvia_u_db=a.lluvia_u,lluvia_d_db=a.lluvia_d,potencia_satelite_w=a.potencia_satelite,diametro_rx_m=a.diametro)
    if a.exportar:
        import matplotlib
        matplotlib.use('Agg')
        r=exportar(e,a.salida)
        fig=figura(e,False); fig.savefig(Path(a.salida)/'panel.png',dpi=160)
        print(json.dumps(r,indent=2,ensure_ascii=False))
    else:
        import matplotlib.pyplot as plt
        figura(e); plt.show()

if __name__=='__main__': main()
