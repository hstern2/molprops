from flask import Flask, request, render_template
from molprops import MolProps

app = Flask(__name__)

h = '''
<!DOCTYPE html>
<head><title>Molecular Properties</title></head>
<body>
<style>
body { font-family: Arial }
table { border-collapse: collapse; }
th, td { border: 1px solid; }
</style>
<p>Enter SMILES, one per line:</p>
<form method="post">
    <textarea name="smiles" rows="10" cols="60"></textarea>
    <br>
    <input type="submit" value="get properties">
    <p></p>
</form>
'''

f = '</body></html>'

@app.route('/', methods=['GET', 'POST'])
def home():
    isPost = request.method == 'POST'
    return h + (MolProps(request.form.get('smiles')) if isPost else '') + f

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=80)
