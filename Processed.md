python3 ../../TuMags_Reduction_Pipeline/process_data_main.py -f YAMLS_pre/config_tumag_01_QSUN_Fe_v0.4.yaml 

- COMM1 
    - DONE Fe (mejorada la v0.4 con cross-interf 8)
    - DONE Mg (mejorada la v0.4 con cross-interf 8)

- 02_EMEF_AR1 
    - DONE Fe (mejorada la v0.4 )
    - DONE Mg (mejorada la v0.4 )

- 02_EMEF_AR2
    - DONE Mg (mejorada la v0.4 )
    - DONE Fe (mejorada la v0.4 )
    - Queda chequear que este todo bien...

- 03_SPOT_SP6
    - Fe DONE
    - Mg DONE (son con código 31, hay que documentarlo y rehacerlo bien)
    
- 04_QS_QS12_1

- 04_QS_QS12_2

- 04_QSUN-QS4ab
    - DONE Fe (mejorada la v0.4 con cross-interf 8)
    - DONE Mg (mejorada la v0.4 con cross-interf 8)
    - Hay datos de Mg que faltan. Habría que añadirlos pero sin repetirlo todo (código 31). Ojo. Sólo si faltan menos de 8!
    
- 04_QSUN-QS8 - started
    - DONE Fe (mejorada la v0.4 con cross-interf 8)
    - DONE Mg (mejorada la v0.4 con cross-interf 8)
    - Hay datos de Mg que faltan. Habría que añadirlos pero sin repetirlo todo (código 31). Ojo. Sólo si faltan menos de 8!

- 05_QSUN-FS1
    - Fe. Hacer una prueba de ROI 4. Hay una network fuerte y puede que se vea afectado por el cross-talk. Necesito chequear los resultados antes. 

- 06_SPOT_AR1 
    - DONE Fe (mejorada la v0.4 )
    - DONE Mg (mejorada la v0.4 )
    - Hay una version _av que tiene el alineamiento por longitud de onda

- 06_SPOT_AR2
PENDIENTE (Faltan datos)

- 08_QSUN
    - QSHC1
    - QSHC2
    - QSHC3







- 11_FLAR_FS1 - DONE
    DONE 1 and 2.06 BUT gradient in Q and last wave in Mg
    Spotted an error in aligment and then in jaegli. ROI regions. 
    repitiendo 0.7 de mg DONE
    repitiendo 1.0 de Mg ...... pintar graficas resultantes DONE
    repitiendo Fe 2.06 DONE
     Procesando Fe PD DONE
     Procesando Mg PD DONE



- 04_QSUN-QS2

- 27_EMEF

- 10_SPOT-SP4a&b

- 10_SPOT-SP9

- 16_LIMB-PL1

- 35_CHOL


- 23_SPOT-SP3&8

02_EMEF, 24_SPOT, y 25_SPOT a la lista de reducción, para Mariarita de parte de Solarnet




./copiar_renombrar_fits.sh /work/obs/TuMag_proccessed_data_2026/11_FLAR_FS1 --only-matching --lv LV_1.0 --from-prefix "11_FLAR_FS1_TM_" --to-prefix "12_FLAR_TM_00_"

