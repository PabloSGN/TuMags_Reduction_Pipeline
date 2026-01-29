Processed:

- 11_FLAR_FS1 - DONE
    DONE 1 and 2.06 BUT gradient in Q and last wave in Mg
    Spotted an error in aligment and then in jaegli. ROI regions. 
    repitiendo 0.7 de mg DONE
    repitiendo 1.0 de Mg ...... pintar graficas resultantes DONE
    repitiendo Fe 2.06 DONE
     Procesando Fe PD DONE
     Procesando Mg PD DONE

- 06_SPOT_AR1
    Calculando Fe 0.5 y 0.7 y 1.0.... 
    Repitiendo 1.0 con jaegli DONE
    PD done en Fe
    - running Mg
    
- COMM1 - DONE
    - procesando hierro sin cont crosst: 
    - Añado _inter8 en label. Ahora tengo que hacer que PD lo vea (modificar código)
    - Voy por el Mg 
    - All done. 

- QS4ab
    - Running config_tumag_04_QSUN_QS4ab_Fe.yaml 


- EMEF



copia comm1 al servidor (toda la carpeta)
copia comm1 a level1
copia comm1 al ordenador

corre 06_SPOT_AR1 de Mg RUNNING

copia spot al servidor (toda la carpeta) y luego a level1


Copia al PC FLARE 1.0 y 1.1 

cp *1.0* /Users/orozco/IdAdA\ Dropbox/David\ orozco\ suárez/Public/datadownload/11_FLAR_FS1 

cp *1.1* /Users/orozco/IdAdA\ Dropbox/David\ orozco\ suárez/Public/datadownload/11_FLAR_FS1 







./copiar_renombrar_fits.sh /work/obs/TuMag_proccessed_data_2026/11_FLAR_FS1 --only-matching --lv LV_1.0 --from-prefix "11_FLAR_FS1_TM_" --to-prefix "12_FLAR_TM_00_"
















