#!/usr/bin/env bash
# Reproduz todos os resultados do manuscrito v5. Executar de qualquer diretório.
# Requisitos: Python>=3.10, numpy, scipy, scikit-learn (ver ../../requirements.txt)
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p results
run(){ out="$1"; shift; echo ">> python3 $*"; python3 "$@" > "results/$out" 2>&1; }

# Galáxias SPARC (Seção 3)
run auditoria_fase1b_log.txt audit2.py
run auditoria_fase1c_log_sigint0.057.txt audit3.py 0.057
# Clusters / BCGs (Seção 1.1)
run clusters_log.txt clusters_test.py
run form_log.txt form_test.py
# Controles de circularidade, ξ_b, previsões (Seções 4-6)
run rgd_v5_log_H67.4.txt rgd_v5.py 67.4
run groups2_log.txt groups_test2.py 67.4
python3 groups_test2.py 70 > results/groups2_log_H70.txt 2>&1
# λ, A em galáxias, monotonicidade (Seções 3, 6, 7)
run rgd_tests_H67.4_log.txt rgd_tests.py 67.4
run rgd_tests_H70_log.txt rgd_tests.py 70
# Exploração anterior (indução radiativa / lenteamento) — opcional
if [ "${1:-}" = "--exploratory" ]; then
  run auditoria_fase1_log.txt audit.py
  run lens_log.txt lens_test.py
  run lens_check_log.txt lens_check.py
  run robust1_A_log.txt robust1.py A
  run robust1_B_log.txt robust1.py B
fi
echo "OK — resultados em v5/code/results/"
