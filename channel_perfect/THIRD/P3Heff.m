 [~, sorted_indices] = sort(h_p_db, 'ascend');
        top_3_indices = sorted_indices(1:3);

        delay_3 = delay(top_3_indices);
        h_p_db_3 = h_p_db(top_3_indices);
        li_3 = li(top_3_indices);
        ki_3 = ki(top_3_indices);
        h_p_3 = h_p(top_3_indices) / sum(h_p(top_3_indices))); % 重新归一化功率
        h_exp_3 = h_exp(top_3_indices);
        P_3 = length(delay_3); % P_3 = 3

        % 2. 计算3径的DD域等效信道 hw3
        hw3 = zeros(M,N);
        for l=0:M-1
            for k=0:N-1
                for i=1:P_3 % P_3 = 3
                    theta=exp(1j*2*pi*ki_3(i)*(li_3(i))/(M*N));
                    delta_term = (l == li_3(i));
                    hw3(l+1,k+1)=hw3(l+1,k+1)+h_p_3(i)*h_exp_3(i)*delta_term*zeta_N(k-ki_3(i),N).*theta;
                end
            end
        end

        % 3. 生成最终的 Heff3
        Heff3 = zeros(M*N, M*N);
        for l = 0:M-1
            for k = 0:N-1
                row_idx = k*M + (l+1);
                for l_prime = 0:M-1
                    for k_prime = 0:N-1
                        col_idx = k_prime*M + (l_prime+1);
                        h_l = mod(l - l_prime, M);
                        h_k = mod(k - k_prime, N);
                        Heff3(row_idx, col_idx) = hw3(h_l + 1, h_k + 1);
                    end
                end
            end
        end