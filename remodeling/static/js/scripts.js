
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
// Função para exibir o modal
function showModal() {
    $('#successModal').modal('show'); // Usa jQuery para mostrar o modal
  }
  
  // Função para lidar com o envio do formulário
  $(document).ready(function() {
    $('form').submit(function(event) {
      event.preventDefault(); // Impede o envio do formulário padrão
  
      var form = $(this);
      var formData = new FormData(form[0]);
  
      $.ajax({
        type: form.attr('method'),
        url: form.attr('action'),
        data: formData,
        cache: false,
        contentType: false,
        processData: false,
        success: function(response) {
          // Se a resposta indicar sucesso, mostra o modal
          if (response.success) {
            showModal();
          }
        }
      });
    });
  });