from flask import Flask,render_template,redirect, url_for,request,redirect,jsonify,Response,send_from_directory,send_file
import os, sys, time,base64
import subprocess
from subprocess import Popen,PIPE
from pymongo import MongoClient
from werkzeug import secure_filename
import numpy
import matplotlib
matplotlib.use('Agg')
from app.module import Expression_profiles_preprocessing
from app.module import Analysis_code

app = Flask(__name__,static_folder='static', static_url_path='')
#app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024i

@app.route('/')
def showRoot():
    return render_template('index.html')

ALLOWED_EXTENSIONS = set(['txt', 'docx'])
def allowed_file(filename):

    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS
def run_function(Case_ID,USER_FOLDER):
    logfile = '>'+USER_FOLDER + 'logfile.txt'
    case = subprocess.Popen(['python3','/var/www/helloworldapp/app/module/Analysis_code.py',Case_ID,USER_FOLDER,logfile], stdout=PIPE, stderr=PIPE, stdin=PIPE) 
    return None

#analysis page

@app.route('/analysis/')
def analysis():
        return render_template('analysis.html')

@app.route('/upload_page', methods=['GET','POST'])
def request_page():
    UPLOAD_FOLDER = '/home/dennistsai/user_upload'
    Case_ID = str(int(time.time()))
    USER_FOLDER = os.path.join(UPLOAD_FOLDER + '/' + Case_ID)
    logfile = USER_FOLDER + 'logfile.txt'
    if (not os.path.exists(USER_FOLDER)): 
       os.mkdir(USER_FOLDER)
    #app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
    upload_files = ['mRNA','lncRNA','miRNA']
         
    for f in upload_files:
        file = request.files[f] 
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            if f == 'mRNA':
               file.save(os.path.join(USER_FOLDER, Case_ID + '_mRNA_file.txt'))
            elif f == 'lncRNA':
               file.save(os.path.join(USER_FOLDER, Case_ID + '_lncRNA_file.txt'))
            else:
               file.save(os.path.join(USER_FOLDER, Case_ID + '_miRNA_file.txt'))
    
    try:
     run_function(Case_ID,USER_FOLDER)
    except Exception as e:
     return str(e)
   
    #case = subprocess.Popen(['python3','/var/www/helloworldapp/app/module/Analysis_code.py',Case_ID,USER_FOLDER,'>',logfile], stdout=PIPE, stderr=PIPE, stdin=PIPE) 
    #Expression_profiles_preprocessing.Concatenating_gene_expression_profile(USER_FOLDER, Case_ID, Case_ID + '_mRNA_file', Case_ID + '_lncRNA_file')    
    
    return url_for('result',case_id=Case_ID)#Case_ID


@app.route('/result/<case_id>')
def result(case_id):
    return render_template('result.html',case_id=case_id)

#check_result page

@app.route('/check_result/')
def check_result():
       return render_template('check_result.html')
    

@app.route('/send_data', methods=['GET','POST'])
def send_data():
    find_ID = request.form['find_ID']
    Cosine_img_path = os.path.join('/home/dennistsai/user_upload',find_ID,find_ID+'_Cosine_lncRNA_association_index_histmatrix_heatmap.png')
    Geometric_img_path = os.path.join('/home/dennistsai/user_upload',find_ID,find_ID+'_Geometric_lncRNA_association_index_histmatrix_heatmap.png')
    Jaccard_img_path = os.path.join('/home/dennistsai/user_upload',find_ID,find_ID+'_Jaccard_lncRNA_association_index_histmatrix_heatmap.png')
    PCC_img_path = os.path.join('/home/dennistsai/user_upload',find_ID,find_ID+'_PCC_lncRNA_association_index_histmatrix_heatmap.png')
    Simpson_img_path = os.path.join('/home/dennistsai/user_upload',find_ID,find_ID+'_Simpson_lncRNA_association_index_histmatrix_heatmap.png')
    
    try:
     subprocess.call(['cp',Cosine_img_path,'/var/www/helloworldapp/app/static/'] )
     subprocess.call(['cp',Geometric_img_path,'/var/www/helloworldapp/app/static/'] )
     subprocess.call(['cp',Jaccard_img_path,'/var/www/helloworldapp/app/static/'] )
     subprocess.call(['cp',PCC_img_path,'/var/www/helloworldapp/app/static/'] )
     subprocess.call(['cp',Simpson_img_path,'/var/www/helloworldapp/app/static/'] )
    except Exception as e:
     return str(e)
    filename = find_ID + '_Cosine_lncRNA_association_index_histmatrix_heatmap.png'
    Cosine_url = url_for('static',filename=filename)
    filename = find_ID + '_Geometric_lncRNA_association_index_histmatrix_heatmap.png'
    Geometric_url = url_for('static',filename=filename)
    filename = find_ID + '_Jaccard_lncRNA_association_index_histmatrix_heatmap.png'
    Jaccard_url = url_for('static',filename=filename)
    filename = find_ID + '_PCC_lncRNA_association_index_histmatrix_heatmap.png'
    PCC_url = url_for('static',filename=filename)
    filename = find_ID + '_Simpson_lncRNA_association_index_histmatrix_heatmap.png'
    Simpson_url = url_for('static',filename=filename)
    
    return jsonify(Cosine_url,Geometric_url,Jaccard_url,PCC_url,Simpson_url,find_ID)
    return url_for('check_result',find_ID=find_ID,image=image)
    
    #return json.dumps({'status':'ok','find_ID':find_ID})
 

if __name__ == '__main__':
    #app.debug = True
    app.run(debug=True())
