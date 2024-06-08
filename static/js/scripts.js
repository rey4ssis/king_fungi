<script>
    // Alterna entre os campos de Texto e Arquivo com base na escolha do usuário
    const opcaoTexto = document.getElementById("opcaoTexto");
    const opcaoArquivo = document.getElementById("opcaoArquivo");
    const campoTexto = document.getElementById("campoTexto");
    const campoArquivo = document.getElementById("campoArquivo");
  
    opcaoTexto.addEventListener("change", () => {
        campoTexto.style.display = "block";
        campoArquivo.style.display = "none";
    });
  
    opcaoArquivo.addEventListener("change", () => {
        campoTexto.style.display = "none";
        campoArquivo.style.display = "block";
    });
  </script>