# beam-fea convenience targets. Override the interpreter with: make run PYTHON=python
PYTHON ?= python3

.PHONY: install test run clean

install:
	$(PYTHON) -m pip install -r requirements.txt

test:
	$(PYTHON) -m pytest -q

run:
	$(PYTHON) scripts/run_all.py

# Removes generated outputs only; keeps results/.gitkeep and ANSYS screenshots.
clean:
	rm -f results/*.csv results/*.png results/summary.md
	rm -rf .pytest_cache
	find . -name "__pycache__" -type d -prune -exec rm -rf {} +
