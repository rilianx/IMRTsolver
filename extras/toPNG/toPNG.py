import png
import json
import sys

SIZE = 20  # Dimension de cada pixel
file_name = "archivo.json"
maximo = 18


def crearMatriz(ancho:int, alto:int):
    """
        Crea una matriz con las dimensiones ingresadas
        como parametro

        Parametros
        ----------
        ancho : int
            Ancho de la matriz que se va a crear
        alto : int
            Alto de la matriz que se va a crear

        Retorna
        -------
        matriz
            Una matriz con las dimensiones especificadas
            como parametro
    """
    
    matriz = []
    for i in range(alto):
        subLista = []
        for j in range(ancho):
            subLista.append(0)
        matriz.append(subLista)
    return matriz

def dibujarRectangulo(matriz, valor, x, y, ancho, alto):
    """
        Dibuja un rectangulo dentro de la matriz ingresada como
        parametro con el valor, posicion y dimension especificada.

        Parametros
        ----------
        x : int
            Posicion en el eje X del rectangulo
        y : int
            Posicion en el eje Y del rectangulo
        valor : int
            Valor que va a tener el rectangulo dentro de la matriz
        ancho : int 
            Ancho del rectangulo que se va a crear
        alto : int
            Alto del rectangulo que se va a crear
    """

    for j in range(y, y+alto):
        for i in range(x, x+ancho):
            matriz[j][i] = valor

def cargarJSON(filename):
    """
        Carga el archivo json con el nombre "filename".
    """
    f = open(filename, "r")
    jString = f.read()
    resultado = json.loads(jString)
    f.close()
    return resultado

def procesarMatriz(matriz, escala, maximo):
    """
        Rescala la matriz y el valor numerico de cada casilla para luego
        generar la imagen.
    """
    matrizR = crearMatriz(len(matriz[0])*escala, len(matriz)*escala)
    escala_color = 255/maximo
    for j in range(len(matriz)):
        for i in range(len(matriz[0])):
            color = max(int(matriz[j][i]*escala_color), 0)
            dibujarRectangulo(matrizR, color, i*escala, j*escala, escala, escala)
    return matrizR

def help():
    print("""
toPNG
-----

    El programa lee un archivo json y genera un png por cada matriz dentro del
    fichero. El formato del archivo es:
    {
        "nombreMatriz" : [[...],[...],[...],...],
        "nombreMatriz2" : [[...],[...],[...],...],
        ...
    }
    Las imagenes se guardaran con el nombre de cada matriz.

Parametros
----------
    -f [file.json]: - Especifica el archivo que se va a leer.
    -m [numero]:    - Especifica el valor maximo que se podra encontrar dentro 
                      de las matrices. Este valor es importante ya que de esto 
                      depende la escala de grices utilizada en la imagen 
                      generada.
    -s [numero]:    - Especifica la escala de la imagen. Por ejemplo, "-s 10"
                      indica que por cada valor dentro de la matriz, se generara
                      un cuadrado de 10x10 pixeles dentro de la imagen.
    """)

if __name__ == "__main__":

    if(len(sys.argv) > 1):
        for i in range(1,len(sys.argv)):
            if(sys.argv[i] == "-f"):
                i+=1
                file_name = sys.argv[i]
            elif(sys.argv[i] == "-m"):
                i+=1
                maximo = int(sys.argv[i])
            elif(sys.argv[i] == "-s"):
                i+=1
                SIZE = int(sys.argv[i])
    else:   
        help()
        exit()
    
    print("Cargando: " + file_name)
    matrices = cargarJSON(file_name)

    # Generando las imagenes
    for nombre, matriz in matrices.items():
        print("Creando: " + nombre + ".png")
        png.from_array(procesarMatriz(matriz, SIZE, maximo), 'L').save(nombre + ".png")

    print("Completado sin problemas")