To PNG
------
El programa lee un archivo json y genera un png por cada matriz dentro del
fichero. El formato del archivo es:
```json
{
        "nombreMatriz" : [[...],[...],[...],...],
        "nombreMatriz2" : [[...],[...],[...],...],
        ...
}
```
Las imagenes se guardaran con el nombre de cada matriz.

Parametros
----------
    -f [file.json]: - Especifica el archivo que se va a leer.
    -m [numero]:    - Especifica el valor maximo que se podría encontrar dentro 
                      de las matrices. Este valor es importante ya que de esto 
                      depende la escala de grices utilizada en la imagen 
                      generada.
    -s [numero]:    - Especifica la escala de la imagen. Por ejemplo, "-s 10"
                      indica que por cada valor dentro de la matriz, se generara
                      un cuadrado de 10x10 pixeles dentro de la imagen.

## Dependencias
- Python3
- Libreria PNG. Para instalar esta librería se puede usar `pip3 install pypng`. 
Mas información en [git](https://github.com/drj11/pypng/).