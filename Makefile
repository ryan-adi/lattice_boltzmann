run:
	python main.py

run_tests:
	python -m unittest discover .

setup: requirements.txt
	pip install -r requirements.txt

clean:
	rm -rf __pycache__