# TuMag header

Description of the information included in the header of the images. 

## Accessing the header. 

Theres is two options to access the header of a specific image:

- 1. Through image ID.

Using the function [read_ID](../image_handler.py#L270) from the module [image_handler](../image_handler.py) you can read the image and obtain the header: 

```python
I, H = image_handler.read_ID("D10-6000")
```

- 2. Through image path.

Using the function [read](../image_handler.py#L90) from the module [image_handler](../image_handler.py) you can read the image and obtain the header: 

```python
I, H = image_handler.read("path/to/image.img")
```

In both cases the header is stored as a dictionary in the variable **H**.

## Header fields.

To access any field use H["field_name"]

| Field name | Description|
|:--------:|:--------:|
| cam | Camera ID. 0 for cam 1, 1 for cam 2. |
| ObservationMode | Observation mode (see the [Observation modes table](../README.md#observation-modes)). |
| nAcc | Number of accumulations. |
| Roix | Size of image in x (default : 2016)  |
| Roiy | Size of image in y (default : 2016) |
| ObservationCounter | Observation counter (see the [Observation Counters](../README. |md#observation-modes))   |
| FW1 | The position of the first filter wheel.  |
| FW2 | The position of the second filter wheel.  |
| hvps_comm_volts | Voltage of the high power supply (etalon) |
| lcvr1_volts | Voltage of first LCVR |
| lcvr2_volts | Volts of second LCVR |
