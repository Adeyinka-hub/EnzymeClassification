#!/usr/bin/python3.6
# Remmy Chen 04/13/2019

import os
import math
import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.naive_bayes import MultinomialNB
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.metrics import classification_report, accuracy_score
from sklearn.metrics import confusion_matrix
pd.options.display.max_columns = 30

PROCESSED_DATA_DIR = 'processed_data'

"""
Loads data into panda frames.
"""
def load_data():
	x_df = pd.read_csv(PROCESSED_DATA_DIR + '/x.csv')
	y_df = pd.read_csv(PROCESSED_DATA_DIR + '/y.csv')
	return x_df, y_df

"""
Separates each A by a space in each AA sequence string.
"""
def preprocess_data(x):
	y = []
	for entry in x:
		y.append(" ".join(entry))
	return y

"""
Extracts data into lists, shuffles data, splits into training 
and test sets, sets cross validation number.
"""
def prep_data(x_df, y_df):
	x = np.array(x_df['AA_seq'].values.tolist())
	x = preprocess_data(x)
	y = np.array(y_df["has_EC_num"].values.tolist())
	assert(len(x) == len(y))
	seed = 42
	split = math.floor(len(x) * 0.8)
	np.random.seed(seed)
	np.random.shuffle(x)
	np.random.seed(seed)
	np.random.shuffle(y)
	x_tr, x_te = np.array(x[:split]), np.array(x[split:])
	y_tr, y_te = np.array(y[:split]), np.array(y[split:])
	cv_num = 5
	return cv_num, x, y, x_tr, x_te, y_tr, y_te




"""
Hyperparameter tuning pipeline.
"""

def tune(cv_num, x_tr, y_tr, x_te):
	pipeline = Pipeline(steps=[
	      ('tfidf', TfidfVectorizer(token_pattern = r"(?u)\b\w+\b")),
	      ('multi_nb', MultinomialNB())
	      ])

	parameters = {'tfidf__ngram_range': [(1, 1), (5, 5), (10, 10), (13, 13), (15, 15), (18, 18), (20, 20)],
	              'multi_nb__alpha': [0, 0.03, 0.07, 0.1, 0.2, 0.5, 1] #(Laplace/Lidstone) smoothing parameter 
	             }

	pipe = GridSearchCV(pipeline, parameters, scoring='accuracy', cv=cv_num, refit='True', return_train_score=True)
	pipe.fit(x_tr, y_tr)
	yhat1_tr = pipe.predict_proba(x_tr)[:, 1]
	yhat1_te = pipe.predict_proba(x_te)[:, 1]
	return pipe, parameters, yhat1_tr, yhat1_te

"""
 Error analysis.
"""
def err_analyze(pipe, parameters, hyperparamType, cv_num):
	fold_train_err_list = []
	fold_test_err_list = []
	for i in range(cv_num):
	    fold_train = []
	    fold_test = []
	    istr = str(i)
	    for j in range(len(parameters[hyperparamType])):
	        f_train = 1 - pipe.cv_results_['split'+istr+'_train_score'][j]
	        f_test = 1 - pipe.cv_results_['split'+istr+'_test_score'][j]
	        fold_train.append(f_train)
	        fold_test.append(f_test)
	    fold_train_err_list.append(fold_train)
	    fold_test_err_list.append(fold_test)

	std_train_list = []
	std_test_list = []
	for i in range(len(parameters[hyperparamType])):
	    std_train_list.append(np.std(np.asarray(fold_train_err_list)[:,i]))
	    std_test_list.append(np.std(np.asarray(fold_test_err_list)[:,i]))
	return fold_train_err_list, fold_test_err_list, std_train_list, std_test_list

if __name__ == "__main__":
	x_df, y_df = load_data()
	cv_num, x, y, x_tr, x_te, y_tr, y_te = prep_data(x_df, y_df)
	##vectorizer = TfidfVectorizer(token_pattern = r"(?u)\b\w+\b") # Tfid filters out strings of length < 2 by default
	##X = vectorizer.fit_transform(x_tr)
	pipe, yhat1_te = tune(cv_num, x_tr, y_tr, x_te)
	print(pipe.best_params_)
	print(pipe.best_score_)
	thresh = 0.5
	accuracy = np.mean(y_te == (yhat1_te >= thresh))
	print(accuracy)








