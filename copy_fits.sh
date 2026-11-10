#!/usr/bin/env bash
# Copia archivos .fits del origen al destino renombrando el prefijo inicial y opcionalmente el nivel LV.
# - Reemplaza el prefijo indicado con --from-prefix por el de --to-prefix SOLO al inicio del nombre.
# - Filtra por LV (p.ej. --lv LV_1.0) si se indica.
# - Permite cambiar la etiqueta del nivel FITS en el archivo resultante mediante --lv_change_to (ej. LV_1.2).
# - Permite especificar el directorio de destino con --dest.
# - No se detiene ante el primer error; sigue y reporta.

set -u  # Fallar si se usan variables no definidas

usage() {
  cat <<'EOF'
Uso:
  copiar_renombrar_fits.sh DIRECTORIOORIGEN \
    --from-prefix "11_FLAR_FS1_TM_" \
    --to-prefix   "12_FLAR_TM_00_" \
    [--dest DIRECTORIODESTINO] \
    [--dry-run] [--overwrite] [--only-matching] \
    [--lv LV_X.Y] [--lv_change_to LV_Z.W]

Obligatorios:
  --from-prefix   Prefijo de entrada que debe estar al principio del nombre.
  --to-prefix     Prefijo de salida que sustituirá al anterior.

Opcionales:
  --dest DIRECTORIO  Directorio de destino donde se copiarán los archivos (Por defecto: level1).
  --lv LV_1.0        Filtra por versión LV a buscar (ej.: LV_1.0, LV_1.1...).
  --lv_change_to LV  Sustituye la etiqueta de nivel encontrada por la indicada (ej.: LV_1.2).
  --dry-run          Muestra qué haría sin copiar.
  --overwrite        Sobrescribe si el archivo destino ya existe.
  --only-matching    Solo procesa archivos cuyo nombre empiece por --from-prefix.
  -h, --help         Muestra esta ayuda.

Notas:
- Renombra prefijo solo si coincide al inicio del nombre de archivo.
- Si se especifica --lv_change_to, sustituye la subcadena de --lv en el nombre de salida.
EOF
}

if [[ $# -lt 1 ]]; then
  usage
  exit 1
fi

ORIG="$1"
shift || true

DRY_RUN=false
OVERWRITE=false
ONLY_MATCHING=false
LV_FILTER=""
LV_CHANGE_TO=""
FROM_PREFIX=""
TO_PREFIX=""
DEST="level1"  # Valor por defecto si no se pasa --dest

# Parseo de argumentos
while (( "$#" )); do
  case "$1" in
    --dry-run) DRY_RUN=true; shift ;;
    --overwrite) OVERWRITE=true; shift ;;
    --only-matching) ONLY_MATCHING=true; shift ;;
    --dest)
      if [[ $# -lt 2 ]]; then
        echo "ERROR: --dest requiere una ruta de directorio" >&2
        exit 1
      fi
      DEST="$2"
      shift 2
      ;;
    --lv)
      if [[ $# -lt 2 ]]; then
        echo "ERROR: --lv requiere un valor (ej.: --lv LV_1.0)" >&2
        exit 1
      fi
      LV_FILTER="$2"
      shift 2
      ;;
    --lv_change_to)
      if [[ $# -lt 2 ]]; then
        echo "ERROR: --lv_change_to requiere un valor (ej.: --lv_change_to LV_1.2)" >&2
        exit 1
      fi
      LV_CHANGE_TO="$2"
      shift 2
      ;;
    --from-prefix)
      if [[ $# -lt 2 ]]; then
        echo "ERROR: --from-prefix requiere un valor (ej.: --from-prefix \"11_FLAR_FS1_TM_\")" >&2
        exit 1
      fi
      FROM_PREFIX="$2"
      shift 2
      ;;
    --to-prefix)
      if [[ $# -lt 2 ]]; then
        echo "ERROR: --to-prefix requiere un valor (ej.: --to-prefix \"12_FLAR_TM_00_\")" >&2
        exit 1
      fi
      TO_PREFIX="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Argumento desconocido: $1" >&2
      usage
      exit 1
      ;;
  esac
done

# Validaciones
if [[ ! -d "$ORIG" ]]; then
  echo "ERROR: El directorio origen no existe o no es un directorio: $ORIG" >&2
  exit 1
fi
if [[ -z "$FROM_PREFIX" || -z "$TO_PREFIX" ]]; then
  echo "ERROR: --from-prefix y --to-prefix son obligatorios." >&2
  usage
  exit 1
fi
if [[ -n "$LV_CHANGE_TO" && -z "$LV_FILTER" ]]; then
  echo "ERROR: --lv_change_to requiere que especifiques también --lv para saber qué sustituir." >&2
  exit 1
fi

mkdir -p "$DEST"

count_processed=0
count_skipped=0
count_conflicts=0
count_errors=0

# Recorremos solo el directorio (sin subdirectorios)
while IFS= read -r -d '' f; do
  base="$(basename "$f")"

  # Filtrado por LV si se pidió
  if [[ -n "$LV_FILTER" ]]; then
    if [[ "$base" != *"_${LV_FILTER}_"* && "$base" != *"${LV_FILTER}_"* && "$base" != *"_${LV_FILTER}"* ]]; then
      ((count_skipped++))
      echo "OMITIDO (no coincide LV '$LV_FILTER'): $base"
      continue
    fi
  fi

  # Renombrado del prefijo solo si coincide al inicio
  if [[ "$base" == ${FROM_PREFIX}* ]]; then
    new="${TO_PREFIX}${base:${#FROM_PREFIX}}"
  else
    if $ONLY_MATCHING; then
      ((count_skipped++))
      echo "OMITIDO (prefijo no coincide): $base"
      continue
    else
      new="$base"  # si no coincide y no pedimos only-matching, procesamos sin cambiar prefijo
    fi
  fi

  # Renombrado de la etiqueta LV en el nombre de salida si se solicitó --lv_change_to
  if [[ -n "$LV_CHANGE_TO" && -n "$LV_FILTER" ]]; then
    new="${new//$LV_FILTER/$LV_CHANGE_TO}"
  fi

  dest_path="$DEST/$new"

  # Conflictos
  if [[ -e "$dest_path" && "$OVERWRITE" == false ]]; then
    ((count_conflicts++))
    echo "CONFLICTO: Ya existe $dest_path (use --overwrite para sobrescribir)."
    continue
  fi

  # Acción
  if $DRY_RUN; then
    echo "DRY-RUN: Copiar '$f' -> '$dest_path'"
  else
    if ! cp -p -- "$f" "$dest_path"; then
      ((count_errors++))
      echo "ERROR: Falló la copia de '$f' a '$dest_path' (continúo con el resto)." >&2
      continue
    fi
    echo "Copiado: $base -> $new"
  fi

  ((count_processed++))
done < <(find "$ORIG" -maxdepth 1 -type f \( -iname '*.fits' -o -iname '*.FITS' \) -print0)

echo
echo "Resumen:"
echo "  Procesados: $count_processed"
echo "  Omitidos:   $count_skipped"
echo "  Conflictos: $count_conflicts"
echo "  Errores:    $count_errors"
echo "  Prefijo FROM: '$FROM_PREFIX'  -> Prefijo TO: '$TO_PREFIX'"
if [[ -n "$LV_CHANGE_TO" ]]; then
  echo "  Nivel FROM:   '$LV_FILTER'  -> Nivel TO:   '$LV_CHANGE_TO'"
fi
echo "  Destino:    $DEST"