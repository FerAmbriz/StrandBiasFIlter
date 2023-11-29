#!/bin/bash
#=========================================#
input_file=$1
#=========================================#

# Prefijo para los archivos de salida
output_prefix="small_file"

# Número de líneas por archivo
lines_per_file=100000

# Leer la primera línea
read -r first_line < "$input_file"

# Contador de archivos
file_counter=1

# Contador de líneas
line_counter=0

while IFS= read -r line
do
  # Si es la primera línea de un nuevo archivo
  if (( line_counter % lines_per_file == 0 ))
  then
    # Si no es el primer archivo, cerrar el archivo anterior
    if (( file_counter > 1 ))
    then
      exec 1>&-
    fi

    # Abrir un nuevo archivo para escritura y escribir la primera línea
    exec 1>"$output_prefix$file_counter.txt"
    echo "$first_line"
  fi

  # Escribir la línea actual en el archivo
  echo "$line"

  # Incrementar contadores
  (( line_counter++ ))
  if (( line_counter % lines_per_file == 0 ))
  then
    (( file_counter++ ))
  fi

done < "$input_file"
