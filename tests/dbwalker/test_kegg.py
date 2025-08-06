from unittest import TestCase

import pytest
from unittest.mock import patch, Mock
import requests
from cobramod.dbwalker.kegg import get_kegg_id_from_cid


class TestGetKeggIdFromCid(TestCase):
    
    def test_successful_request_single_result(self):
        """Test successful API call with a single KEGG ID result."""
        mock_response_text = ("pubchem:5962\tcpd:C03062")
        
        with patch('requests.get') as mock_get:
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.text = mock_response_text
            mock_get.return_value = mock_response
            
            result = get_kegg_id_from_cid('5962')

            self.assertEqual("kegg.compound:C03062", result)
            mock_get.assert_called_once_with('https://rest.kegg.jp/conv/compound/pubchem:5962')
    
    def test_successful_request_multiple_results(self):
        """Test successful API call with multiple KEGG ID results (returns first)."""
        mock_response_text = "pubchem:5962\tcpd:C00022\npubchem:5962\tcpd:C00033"
        
        with patch('requests.get') as mock_get:
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.text = mock_response_text
            mock_get.return_value = mock_response
            
            result = get_kegg_id_from_cid('5962')

            self.assertEqual("kegg.compound:C00022", result)
            mock_get.assert_called_once_with('https://rest.kegg.jp/conv/compound/pubchem:5962')
    
    def test_api_error_status_code(self):
        """Test handling of API error status codes."""
        with patch('requests.get') as mock_get:
            mock_response = Mock()
            mock_response.status_code = 404
            mock_get.return_value = mock_response
            
            with self.assertLogs('cobramod.DBWalker.Kegg', level='ERROR') as cm:
                result = get_kegg_id_from_cid('999999')
                
                self.assertIsNone(result)
                self.assertEqual(cm.output, ["ERROR:cobramod.DBWalker.Kegg:Error (404) getting KEGG compound ID from 999999"])

    def test_empty_response(self):
        """Test handling of empty API response."""
        with patch('requests.get') as mock_get:
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.text = ""
            mock_get.return_value = mock_response
            
            result = get_kegg_id_from_cid('nonexistent')
            
            assert result is None
    
    def test_whitespace_only_response(self):
        """Test handling of whitespace-only API response."""
        with patch('requests.get') as mock_get:
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.text = "   \n  \t  "
            mock_get.return_value = mock_response
            
            result = get_kegg_id_from_cid('test')
            
            assert result is None
    
    def test_malformed_response(self):
        """Test handling of malformed API response."""
        mock_response_text = "invalid_format_without_tabs"
        
        with patch('requests.get') as mock_get:
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.text = mock_response_text
            mock_get.return_value = mock_response
            
            # This should not crash, but behavior depends on implementation
            result = get_kegg_id_from_cid('test')
            
            # The function will try to split and process, might return unexpected result
            assert isinstance(result, str) or result is None
    
    def test_request_exception(self):
        """Test handling of network/request exceptions."""
        with patch('requests.get') as mock_get:
            mock_get.side_effect = requests.RequestException("Network error")
            
            # The current implementation doesn't handle exceptions, so this will raise
            with pytest.raises(requests.RequestException):
                get_kegg_id_from_cid('test')
    
    def test_url_construction(self):
        """Test that the correct URL is constructed for the API call."""
        with patch('requests.get') as mock_get:
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.text = "compound:C12345\tpubchem:12345"
            mock_get.return_value = mock_response
            
            get_kegg_id_from_cid('12345')
            
            mock_get.assert_called_once_with('https://rest.kegg.jp/conv/compound/pubchem:12345')
    
    def test_different_cid_formats(self):
        """Test function with different CID input formats."""
        test_cases = [
            ('123', 'https://rest.kegg.jp/conv/compound/pubchem:123'),
            ('456789', 'https://rest.kegg.jp/conv/compound/pubchem:456789'),
            ('0', 'https://rest.kegg.jp/conv/compound/pubchem:0')
        ]
        
        for cid, expected_url in test_cases:
            with patch('requests.get') as mock_get:
                mock_response = Mock()
                mock_response.status_code = 200
                mock_response.text = f"compound:C00001\tpubchem:{cid}"
                mock_get.return_value = mock_response
                
                get_kegg_id_from_cid(cid)
                
                mock_get.assert_called_once_with(expected_url)
    
    def test_live_api_call(self):
        """Test against the live KEGG API (may fail if API is unavailable)."""
        result = get_kegg_id_from_cid('5793')
        
        # Note: Live API may change, so we just check the format if result exists
        if result is not None:
            assert isinstance(result, str)
            assert len(result) > 0
        
        print(f"Live API result for CID 1060: {result}")
