Processed:

- 11_FLAR_FS1 - DONE
    DONE 1 and 2.06 BUT gradient in Q and last wave in Mg
    Spotted an error in aligment and then in jaegli. ROI regions. 
    repitiendo 0.7 de mg DONE
    repitiendo 1.0 de Mg ...... pintar graficas resultantes DONE
    repitiendo Fe 2.06 DONE
     Procesando Fe PD DONE
     Procesando Mg PD DONE

- 06_SPOT_AR1 - DONE
    Calculando Fe 0.5 y 0.7 y 1.0.... 
    Repitiendo 1.0 con jaegli DONE
    PD done en Fe
    - running Mg

- COMM1 - DONE
    - procesando hierro sin cont crosst: 
    - Añado _inter8 en label. Ahora tengo que hacer que PD lo vea (modificar código)
    - Voy por el Mg 
    - All done. 

- 04_QSUN-QS4ab
    - Running config_tumag_04_QSUN_QS4ab_Fe.yaml 
    - Running config_tumag_04_QSUN_QS4ab_Mg.yaml 
    Falta los ninter_8 de los dos.

2026-01-30 20:26:27,830 [SpawnProcess-1] INFO:  processing ocs: 18 ........... 
2026-01-30 20:26:27,830 [SpawnProcess-1] INFO:  processing ocs: 18 which corresponds to 1 with 17 total images
2026-01-30 20:26:27,830 [SpawnProcess-1] ERROR:   >> Error. The ocs: 18 number of images 17 does not coincide with the timeline info: 80
2026-01-30 20:26:51,223 [SpawnProcess-6] INFO:  processing ocs: 19 ........... 
2026-01-30 20:26:51,223 [SpawnProcess-6] INFO:  processing ocs: 19 which corresponds to 1 with 80 total images


- 04_QSUN-QS8

- 04_QSUN-QS2

- 27_EMEF

- 02_EMEF

- 10_SPOT-SP4a&b

- 10_SPOT-SP9

- 16_LIMB-PL1

- 35_CHOL

- 05_QSUN-FS_1

- 23_SPOT-SP3&8



Copia al PC FLARE 1.0 y 1.1 

cp *1.0* /Users/orozco/IdAdA\ Dropbox/David\ orozco\ suárez/Public/datadownload/11_FLAR_FS1 

cp *1.1* /Users/orozco/IdAdA\ Dropbox/David\ orozco\ suárez/Public/datadownload/11_FLAR_FS1 







./copiar_renombrar_fits.sh /work/obs/TuMag_proccessed_data_2026/11_FLAR_FS1 --only-matching --lv LV_1.0 --from-prefix "11_FLAR_FS1_TM_" --to-prefix "12_FLAR_TM_00_"

