# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 17:25:53 2017

@author: User
"""

import requests
def test():
	user_info = {'name': 'letian', 'password': '123'}
	r = requests.post("secure-garden-72131.herokuapp.com/register", data=user_info)

	return r.text
