import os
import numpy as np
import pydicom
from scipy.interpolate import griddata
from utils import logger

class BaseFileHandler:
    """다양한 파일 핸들러를 위한 기본 클래스"""
    def __init__(self):
        self.filename = None
        self.pixel_data = None
        self.phys_x_mesh = None
        self.phys_y_mesh = None
        self.physical_extent = None
        self.origin_x = 0
        self.origin_y = 0
        self.pixel_spacing = 1.0
        self.dose_bounds = None
        self.crop_pixel_offset = (0, 0)
        
    def get_filename(self):
        """파일명 반환"""
        if self.filename:
            return os.path.basename(self.filename)
        return None
    
    def get_physical_extent(self):
        """물리적 좌표 범위 반환"""
        return self.physical_extent
    
    def get_origin_coords(self):
        """원점 좌표 반환"""
        return self.origin_x, self.origin_y
        
    def get_spacing(self):
        """픽셀 간격 반환"""
        return self.pixel_spacing, self.pixel_spacing
    
    def get_pixel_data(self):
        """픽셀 데이터 반환"""
        return self.pixel_data
    
    def create_physical_coordinates(self):
        """물리적 좌표계 생성 (추상 메서드)"""
        raise NotImplementedError("서브클래스에서 구현해야 함")
    
    def physical_to_pixel_coord(self, phys_x, phys_y):
        """물리적 좌표(mm)를 픽셀 좌표로 변환 (추상 메서드)"""
        raise NotImplementedError("서브클래스에서 구현해야 함")
    
    def pixel_to_physical_coord(self, pixel_x, pixel_y):
        """픽셀 좌표를 물리적 좌표(mm)로 변환 (추상 메서드)"""
        raise NotImplementedError("서브클래스에서 구현해야 함")

    def _calculate_bounds_from_mask(self, mask, margin_mm=0):
        if not np.any(mask):
            return None

        rows = np.any(mask, axis=1)
        cols = np.any(mask, axis=0)

        row_indices = np.where(rows)[0]
        col_indices = np.where(cols)[0]

        min_row, max_row = row_indices[0], row_indices[-1]
        min_col, max_col = col_indices[0], col_indices[-1]

        min_phys_x, min_phys_y = self.pixel_to_physical_coord(min_col, min_row)
        max_phys_x, max_phys_y = self.pixel_to_physical_coord(max_col, max_row)

        if margin_mm > 0:
            min_phys_x -= margin_mm
            max_phys_x += margin_mm
            min_phys_y -= margin_mm
            max_phys_y += margin_mm

        bounds = {
            'min_x': min_phys_x, 'max_x': max_phys_x,
            'min_y': min_phys_y, 'max_y': max_phys_y
        }
        logger.info(f"계산된 선량 경계 (물리적 좌표, {margin_mm}mm 여백 포함): {bounds}")
        return bounds
        
    def open_file(self, filename):
        """파일 로드 (추상 메서드)"""
        try:
            self.filename = filename
            # 파일 로드 로직 (서브클래스에서 구현)
            return True
        except Exception as e:
            error_msg = f"파일 로드 오류: {str(e)}"
            logger.error(error_msg)
            return False, error_msg


class DicomFileHandler(BaseFileHandler):
    """DICOM RT dose 파일을 처리하는 클래스"""
    def __init__(self):
        super().__init__()
        self.dicom_data = None
        self.dicom_origin_x = 0
        self.dicom_origin_y = 0
        self.pixel_spacing = 1.0
        
    def open_file(self, filename):
        """DICOM RT dose 파일 로드"""
        try:
            self.filename = filename
            self.dicom_data = pydicom.dcmread(filename)
            
            if self.dicom_data.Modality != 'RTDOSE':
                raise ValueError("Selected file is not an RT Dose file")
                
            self.pixel_data = self.dicom_data.pixel_array * self.dicom_data.DoseGridScaling
            
            if hasattr(self.dicom_data, 'PixelSpacing'):
                self.pixel_spacing = float(self.dicom_data.PixelSpacing[0])
            
            height, width = self.pixel_data.shape
            
            if hasattr(self.dicom_data, 'ImagePositionPatient'):
                # ImagePositionPatient: [x, z, y]
                # 물리적 위치를 픽셀 간격으로 나누어 픽셀 원점 계산
                self.dicom_origin_x = int(round(self.dicom_data.ImagePositionPatient[0] / self.pixel_spacing))
                self.dicom_origin_y = int(round(self.dicom_data.ImagePositionPatient[2] / self.pixel_spacing)) # y좌표는 index 2
                
                logger.info(f"DICOM 픽셀 원점 설정 (from ImagePositionPatient): x={self.dicom_origin_x}, y={self.dicom_origin_y}")
                logger.info(f"사용된 ImagePositionPatient 값: x={self.dicom_data.ImagePositionPatient[0]}, y={self.dicom_data.ImagePositionPatient[2]}")
            else:
                self.dicom_origin_x = -width // 2
                self.dicom_origin_y = -height // 2
                logger.warning(f"DICOM 원점 정보 없음. 이미지 중심으로 기본값 설정: ({self.dicom_origin_x}, {self.dicom_origin_y})")
                
            self.create_physical_coordinates()

            # 원본 데이터 및 좌표 보존
            full_pixel_data = self.pixel_data
            full_phys_x_mesh = self.phys_x_mesh
            full_phys_y_mesh = self.phys_y_mesh

            self.dose_bounds = self.calculate_dose_bounds(margin_mm=20) # 2cm 여백으로 ROI 자동 크롭

            if self.dose_bounds:
                bounds = self.dose_bounds
                min_px = int(round(bounds['min_x'] / self.pixel_spacing - self.dicom_origin_x))
                max_px = int(round(bounds['max_x'] / self.pixel_spacing - self.dicom_origin_x))
                min_py = int(round(bounds['min_y'] / self.pixel_spacing - self.dicom_origin_y))
                max_py = int(round(bounds['max_y'] / self.pixel_spacing - self.dicom_origin_y))

                if min_py > max_py: min_py, max_py = max_py, min_py
                if min_px > max_px: min_px, max_px = max_px, min_px

                h, w = full_pixel_data.shape
                min_py, max_py = max(0, min_py), min(h, max_py + 1)
                min_px, max_px = max(0, min_px), min(w, max_px + 1)

                self.pixel_data = full_pixel_data[min_py:max_py, min_px:max_px]
                self.phys_x_mesh = full_phys_x_mesh[min_py:max_py, min_px:max_px]
                self.phys_y_mesh = full_phys_y_mesh[min_py:max_py, min_px:max_px]
                self.crop_pixel_offset = (min_px, min_py)
                self.physical_extent = [self.phys_x_mesh.min(), self.phys_x_mesh.max(), self.phys_y_mesh.min(), self.phys_y_mesh.max()]
                logger.info(f"DICOM data has been cropped to ROI. New shape: {self.pixel_data.shape}")
            
            return True
            
        except Exception as e:
            error_msg = f"DICOM 파일 로드 오류: {e}"
            logger.error(error_msg)
            return False, error_msg

    def calculate_dose_bounds(self, dicom_image=None, threshold_percent=0, margin_mm=0):
        """선량 임계값 또는 0이 아닌 선량 영역에 대한 경계를 계산합니다."""
        if dicom_image is None:
            dicom_image = self.get_pixel_data()
        if dicom_image is None:
            return None

        if threshold_percent > 0:
            max_dose = np.max(dicom_image)
            threshold_val = (threshold_percent / 100.0) * max_dose
            mask = dicom_image >= threshold_val
        else:
            mask = dicom_image > 0
        
        return self._calculate_bounds_from_mask(mask, margin_mm)

    def create_physical_coordinates(self):
        if self.pixel_data is None: return
        height, width = self.pixel_data.shape
        phys_x = (np.arange(width) + self.dicom_origin_x) * self.pixel_spacing
        phys_y = (np.arange(height) + self.dicom_origin_y) * self.pixel_spacing
        self.phys_x_mesh, self.phys_y_mesh = np.meshgrid(phys_x, phys_y)
        self.physical_extent = [phys_x.min(), phys_x.max(), phys_y.min(), phys_y.max()]

    def physical_to_pixel_coord(self, phys_x, phys_y):
        full_grid_px = phys_x / self.pixel_spacing - self.dicom_origin_x
        full_grid_py = phys_y / self.pixel_spacing - self.dicom_origin_y
        
        cropped_px = int(round(full_grid_px - self.crop_pixel_offset[0]))
        cropped_py = int(round(full_grid_py - self.crop_pixel_offset[1]))
        
        return cropped_px, cropped_py
    
    def pixel_to_physical_coord(self, pixel_x, pixel_y):
        full_grid_px = pixel_x + self.crop_pixel_offset[0]
        full_grid_py = pixel_y + self.crop_pixel_offset[1]
        
        phys_x = (full_grid_px + self.dicom_origin_x) * self.pixel_spacing
        phys_y = (full_grid_py + self.dicom_origin_y) * self.pixel_spacing
        return phys_x, phys_y

    def get_origin_coords(self):
        return self.dicom_origin_x, self.dicom_origin_y

    def get_spacing(self):
        return self.pixel_spacing, self.pixel_spacing

    def get_patient_info(self):
        """환자 정보(기관명, ID, 이름) 반환"""
        if self.dicom_data is None:
            return None, None, None

        institution = self.dicom_data.get("InstitutionName", "N/A")
        patient_id = self.dicom_data.get("PatientID", "N/A")
        patient_name = self.dicom_data.get("PatientName", "N/A")

        return institution, patient_id, patient_name

class MCCFileHandler(BaseFileHandler):
    """MCC 파일을 처리하는 클래스"""
    def __init__(self):
        super().__init__()
        self.matrix_data = None
        self.device_type = None
        self.task_type = None
        self.n_rows = None
        self.mcc_origin_x = 0
        self.mcc_origin_y = 0
        self.mcc_spacing_x = 1.0
        self.mcc_spacing_y = 1.0
                
    def get_matrix_data(self):
        return self.matrix_data

    def get_interpolated_matrix_data(self, method='cubic'):
        if self.matrix_data is None: return None
        data = self.matrix_data.copy()
        data[data < 0] = np.nan
        valid_points_indices = np.where(~np.isnan(data))
        if len(valid_points_indices[0]) < 4: return self.matrix_data 
        grid_y, grid_x = np.mgrid[0:data.shape[0], 0:data.shape[1]]
        interpolated_data = griddata(
            np.array(list(zip(valid_points_indices[0], valid_points_indices[1]))),
            data[valid_points_indices],
            (grid_y, grid_x),
            method=method
        )
        return interpolated_data        
    
    def open_file(self, filename):
        try:
            self.filename = filename
            with open(filename, "r") as file:
                lines = file.read().split()
            
            self.device_type, self.task_type = self.detect_device_type(lines)
            
            N_begin = lines.count("BEGIN_DATA")
            self.n_rows = N_begin
            self.matrix_data = self.extract_data(lines, N_begin, self.device_type, self.task_type)
            self.pixel_data = self.matrix_data
            
            self._set_device_parameters()
            self.create_physical_coordinates()
            
            logger.info(f"MCC 파일 로드 완료: {self.get_device_name()}")
            return True
                        
        except Exception as e:
            error_msg = f"File open error: {str(e)}"
            logger.error(error_msg)
            return False, error_msg

    def crop_to_bounds(self, bounds):
        """Crops the MCC data based on physical coordinate bounds from DICOM."""
        if self.matrix_data is None or bounds is None:
            return

        full_matrix_data = self.matrix_data
        full_phys_x_mesh = self.phys_x_mesh
        full_phys_y_mesh = self.phys_y_mesh

        # MCC y-axis is inverted. min_phys_y maps to max_pixel_y.
        min_px, max_py_from_min_y = self.physical_to_pixel_coord(bounds['min_x'], bounds['min_y'])
        max_px, min_py_from_max_y = self.physical_to_pixel_coord(bounds['max_x'], bounds['max_y'])

        min_py, max_py = min_py_from_max_y, max_py_from_min_y

        if min_py > max_py: min_py, max_py = max_py, min_py
        if min_px > max_px: min_px, max_px = max_px, min_px

        h, w = full_matrix_data.shape
        min_py, max_py = max(0, min_py), min(h, max_py + 1)
        min_px, max_px = max(0, min_px), min(w, max_px + 1)

        self.matrix_data = full_matrix_data[min_py:max_py, min_px:max_px]
        self.pixel_data = self.matrix_data
        self.phys_x_mesh = full_phys_x_mesh[min_py:max_py, min_px:max_px]
        self.phys_y_mesh = full_phys_y_mesh[min_py:max_py, min_px:max_px]
        self.crop_pixel_offset = (min_px, min_py)
        self.physical_extent = [self.phys_x_mesh.min(), self.phys_x_mesh.max(), self.phys_y_mesh.min(), self.phys_y_mesh.max()]
        
        logger.info(f"MCC data has been cropped to DICOM ROI. New shape: {self.matrix_data.shape}")

    def calculate_dose_bounds(self, mcc_image=None, margin_mm=0):
        """0이 아닌 선량 영역에 대한 경계를 계산합니다."""
        if mcc_image is None:
            mcc_image = self.get_pixel_data()
        if mcc_image is None:
            return None

        mask = mcc_image > 0
        return self._calculate_bounds_from_mask(mask, margin_mm)
    
    def _set_device_parameters(self):
        if self.device_type == 2:  # 1500
            self.mcc_origin_x = 27
            self.mcc_origin_y = 27
            self.mcc_spacing_x = 5.0
            self.mcc_spacing_y = 5.0
        else:  # 725
            self.mcc_origin_x = 13
            self.mcc_origin_y = 13
            self.mcc_spacing_x = 10.0
            self.mcc_spacing_y = 10.0
    
    def extract_data(self, lines, N_begin, device_type, task_type):
        # ... (이전과 동일한 데이터 추출 로직)
        try:
            if task_type == 1: oct_read_intv = 3
            else: oct_read_intv = 2
            st_idx_j, end_idx_j = [], []
            for j, line in enumerate(lines):
                if line == "BEGIN_DATA": st_idx_j.append(j)
                elif line == "END_DATA": end_idx_j.append(j)
            delimt_ind = np.zeros((len(st_idx_j), 3), dtype=int)
            delimt_ind[:, 0], delimt_ind[:, 1] = st_idx_j, end_idx_j
            delimt_ind[:, 2] = delimt_ind[:, 1] - delimt_ind[:, 0] - 1
            matrix_octavius_mat_tmp = np.zeros((N_begin, N_begin)) - 1
            x_lngt = int(delimt_ind[0, 2] / oct_read_intv)
            if device_type == 2 and task_type == 1:
                for j in range(0, N_begin, 2):
                    for k in range(x_lngt-1):
                        matrix_octavius_mat_tmp[j, 2*k+1] = float(lines[delimt_ind[j, 0] + 2 + oct_read_intv*(k)])
                for j in range(1, N_begin, 2):
                    for k in range(1, x_lngt):
                        matrix_octavius_mat_tmp[j, 2*k] = float(lines[delimt_ind[j, 0] + 2 + oct_read_intv*(k-1)])
            else:
                for j in range(N_begin):
                    for k in range(x_lngt):
                        matrix_octavius_mat_tmp[j, k] = float(lines[delimt_ind[j, 0] + 2 + oct_read_intv*k])
            return matrix_octavius_mat_tmp
            # return np.flipud(matrix_octavius_mat_tmp)
        except Exception as e:
            logger.error(f"Data extraction error: {str(e)}")
            raise
        
    def detect_device_type(self, lines):
        try:
            is_1500 = "SCAN_DEVICE=OCTAVIUS_1500_XDR" in lines
            is_merged = "SCAN_OFFAXIS_CROSSPLANE=0.00" in lines
            return (2 if is_1500 else 1), (2 if is_merged else 1)
        except Exception as e:
            logger.error(f"Device type detection error: {str(e)}")
            raise

    def get_device_name(self):
        if self.device_type == 2: return "OCTAVIUS 1500" + (" with merge" if self.task_type == 2 else "")
        else: return "OCTAVIUS 725" + (" with merge" if self.task_type == 2 else "")
                
    def get_origin_coords(self):
        return self.mcc_origin_x, self.mcc_origin_y
            
    def get_spacing(self):
        return self.mcc_spacing_x, self.mcc_spacing_y

    def create_physical_coordinates(self):
        if self.matrix_data is None: return
        height, width = self.matrix_data.shape
        phys_x = (np.arange(width) - self.mcc_origin_x) * self.mcc_spacing_x
        phys_y = -(np.arange(height) - self.mcc_origin_y) * self.mcc_spacing_y  # y축 반전
        self.phys_x_mesh, self.phys_y_mesh = np.meshgrid(phys_x, phys_y)
        self.physical_extent = [phys_x.min(), phys_x.max(), phys_y.min(), phys_y.max()]
            
    def physical_to_pixel_coord(self, phys_x, phys_y):
        full_grid_px = phys_x / self.mcc_spacing_x + self.mcc_origin_x
        full_grid_py = -(phys_y / self.mcc_spacing_y) + self.mcc_origin_y
        
        cropped_px = int(round(full_grid_px - self.crop_pixel_offset[0]))
        cropped_py = int(round(full_grid_py - self.crop_pixel_offset[1]))
        
        return cropped_px, cropped_py
    
    def pixel_to_physical_coord(self, pixel_x, pixel_y):
        full_grid_px = pixel_x + self.crop_pixel_offset[0]
        full_grid_py = pixel_y + self.crop_pixel_offset[1]
        
        phys_x = (full_grid_px - self.mcc_origin_x) * self.mcc_spacing_x
        phys_y = -(full_grid_py - self.mcc_origin_y) * self.mcc_spacing_y
        return phys_x, phys_y
