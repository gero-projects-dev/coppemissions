
import streamlit as st
import subprocess
import sys
from pathlib import Path

st.set_page_config(
    page_title="COPPEmissions",
    page_icon="🌱",
    layout="wide",
)

st.title("COPPEmissions")
st.subheader("Estimativa de emissões de gases de efeito estufa")

st.write(
    "Esta aplicação permite executar os casos de estudo "
    "disponíveis no modelo COPPEmissions."
)

PROJECT_DIR = Path(__file__).resolve().parent

casos = {
    "Expansão e reinjeção": (
        "cases/main_buzios_expreinj_novo_final.py"
    ),
    "Reinjeção total": (
        "cases/main_buzios_reinjtotal_novo_final.py"
    ),
}

caso = st.selectbox(
    "Selecione o caso de estudo:",
    list(casos.keys()),
)

if st.button("Executar modelo", type="primary"):
    arquivo = PROJECT_DIR / casos[caso]

    if not arquivo.exists():
        st.error(f"Arquivo não encontrado: {arquivo}")
    else:
        st.info("Executando o modelo...")

        resultado = subprocess.run(
            [sys.executable, str(arquivo)],
            cwd=PROJECT_DIR,
            capture_output=True,
            text=True,
            timeout=600,
        )

        if resultado.returncode == 0:
            st.success("Execução concluída!")
            if resultado.stdout:
                st.text_area(
                    "Saída do modelo",
                    resultado.stdout,
                    height=300,
                )
        else:
            st.error("O modelo encontrou um erro.")
            st.code(resultado.stderr or resultado.stdout)
