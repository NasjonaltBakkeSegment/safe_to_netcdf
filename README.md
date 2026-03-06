# Data Transformation Tool

This repository contains a Python script designed to transform Copernicus Sentinel products in **Standard Archive Format for Europe (SAFE)** into either **netCDF** or **GeoTIFF** formats. The script supports Sentinel-1 GRD, Sentinel-2, and Sentinel-3 OLCI data.

> Note: transformation to **GeoTIFF** is not currently available

---

## Features

- **Supports Sentinel Missions**:
  - Sentinel-1 (`S1`) - GRD data only
  - Sentinel-2 (`S2`)
  - Sentinel-3 (`S3`) - OLCI data only
- **Output Formats**:
  - `netCDF`
  - `GeoTIFF`
- **Key Capabilities**:
  - Handles SAFE files for Sentinel missions.
  - Generates and updates global and variable attributes for scientific datasets.
  - Adds specialized metadata and calibration layers to the output files.
  - Provides logging and memory profiling during execution.

---

## Usage

Run the script using:

```bash
python main.py <product_name> [options]
```

## Example

```bash
PRODUCT_NAME="S2A_MSIL2A_20260216T142251_N0512_R096_T25WET_20260216T202508"
SAFEDIR="/path/to/safe_directory"  # SAFE directory
NETCDFDIR="/path/to/netcdf_directory"  # netCDF directory
GEOTIFFDIR="/path/to/geotiff_directory"  # GeoTIFF directory
TMPDIR="/path/to/temporary_directory"  # Temporary directory
GLOBAL_ATTRIBUTES_CONFIG="/path/to/config/global_attributes.yaml"  # Global attributes config
VARIABLE_ATTRIBUTES_CONFIG="/path/to/config/variable_attributes.yaml"  # Variable attributes config
PRODUCT_ID="59d354b3-8298-4675-beef-72b756ebb3e5"  # Product ID
TARGET="netCDF"  # Target format, either 'netCDF' or 'geoTIFF'

# Run the Python script with the defined parameters
python3 transform.py \
    "$PRODUCT_NAME" \
    --safedir "$SAFEDIR" \
    --netcdfdir "$NETCDFDIR" \
    --geotiffdir "$GEOTIFFDIR" \
    --tmpdir "$TMPDIR" \
    --global_attributes_config "$GLOBAL_ATTRIBUTES_CONFIG" \
    --variable_attributes_config "$VARIABLE_ATTRIBUTES_CONFIG" \
    --product_id "$PRODUCT_ID" \
    --target "$TARGET"
```

## Input Parameters

The script takes the following arguments:

### Positional Argument:
- `product_name` *(str)*: Name of the product to be transformed, without extension.

### Optional Arguments:
- `--safedir` *(str)*: Path to the directory containing the SAFE files.
- `--netcdfdir` *(str)*: Path to the output directory for netCDF files - or input if netCDF to geotiff.
- `--geotiffdir` *(str)*: Path to the output directory for GeoTIFF files.
- `--tmpdir` *(str)* (default: `tmp/`): Path where to save the uncompressed zip data temporarily. Deleted after usage.
- `--global_attributes_config` *(str)*: Path to the global attributes configuration YAML file (default: `config/global_attributes.yaml`).
- `--variable_attributes_config` *(str)*: Path to the variable attributes configuration YAML file (default: `config/variable_attributes.yaml`).
- `--product_id` *(str)*: Optional product ID to include in the global attributes.
- `--target` *(str)*: Output format, either `netCDF` or `GeoTIFF` (default: `netCDF`).