from flask import Flask, request, render_template
from molprops import MolProps, MolPlots

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

hp = '''
<!DOCTYPE html>
<head><title>Molecular Properties - Plots</title></head>
<body>
<style>
body { font-family: Arial }
table { border-collapse: collapse; }
th, td { border: 1px solid; }
</style>
<p>Enter group labels beginning with #, followed by SMILES, one per line. For example,</p>
<pre>
# straight-chain alkanes
CCC
CCCCC

# alcohols
CCO
CC(C)CO
</pre>
<form method="post">
    <textarea name="smiles" rows="10" cols="60"></textarea>
    <br>
    <input type="submit" value="get properties">
    <p></p>
</form>
'''

f = '</body></html>'

def isPost():
    return request.method == 'POST'

@app.route('/', methods=['GET', 'POST'])
def home():
    return h + (MolProps(request.form.get('smiles')) if isPost() else '') + f

@app.route('/plots', methods=['GET', 'POST'])
def plots():
    return hp + (MolPlots(request.form.get('smiles')) if isPost() else '') + f

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=80)
