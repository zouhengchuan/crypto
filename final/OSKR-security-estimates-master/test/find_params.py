from OSKR import OSKR_to_MLWE, OSKRParameterSet
from OSKR_failure import error_probability
from MLWE_security import MLWE_summarize_attacks

def find_parameters(min_security_strength, min_error_rate_exp):
    # 定义参数搜索范围
    eta_k_range = range(2, 5)  # eta_k 的可能取值范围
    eta_e_range = range(2, 5)  # eta_e 的可能取值范围
    u_range = range(12, 8, -1)  # u 的可能取值范围
    v_range = range(7, 2, -1)  # v 的可能取值范围
    
    # 存储满足条件的参数组合
    valid_params = []

    # 搜索参数空间
    for eta_k in eta_k_range:
        for eta_e in eta_e_range:
            for eta_e_ct in range(eta_e, 5):
                for u in u_range:
                    for v in v_range:
                        # 创建参数集
                        ps_test = OSKRParameterSet(384, 2, eta_k, eta_e, 3457, 2**u, 2**v, ke_ct = eta_e_ct)
                        
                        # 计算安全强度
                        security_strength = MLWE_summarize_attacks(OSKR_to_MLWE(ps_test))
                        
                        # 计算错误率
                        error_rate = error_probability(ps_test)
                        
                        # 检查是否满足条件
                        if security_strength > min_security_strength and error_rate < -min_error_rate_exp:
                            valid_params.append((eta_k, eta_e, eta_e_ct, u, v))
    
    return valid_params

# 调用函数并传入所需的最小安全强度和错误率指数
valid_params = find_parameters(160, 125)
print("Valid parameters:", valid_params)