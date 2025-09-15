"""
This module provides utility functions for the application, including logging and array manipulation.
"""
import os
import logging
import numpy as np
from datetime import datetime

# Check and create log directory
log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logs')
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

# 로그 파일 설정 (일별 로그 파일)
log_file = os.path.join(log_dir, f'gamma_analysis_{datetime.now().strftime("%Y%m%d")}.log')

# 로거 설정
def setup_logger(name):
    """Sets up a logger for the application.

    Args:
        name (str): The name of the logger.

    Returns:
        logging.Logger: The configured logger object.
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    # 이미 핸들러가 추가되어 있다면 추가하지 않음
    if not logger.handlers:
        # 파일 핸들러
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        file_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(file_format)
        logger.addHandler(file_handler)
        
        # 콘솔 핸들러
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(file_format)
        logger.addHandler(console_handler)
    
    return logger

# 기본 로거 생성
logger = setup_logger('gamma_analysis')

def find_nearest_index(array, value):
    """Returns the index of the element in the array that is closest to the given value.

    Args:
        array (numpy.ndarray): The array to search.
        value (float): The value to find.

    Returns:
        int: The index of the nearest value.
    """
    return np.argmin(np.abs(array - value))
