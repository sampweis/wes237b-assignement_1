#include "huffman.h"

#include <math.h>
#include <string.h>
#include <limits.h>
#include <queue>

class huffman_codec
{

private:
	// huffman tree node
	struct node
	{
		// for leaf nodes, contains the character indicated by the node, null otherwise
		unsigned char character;

		// for leaf nodes, the frequency of the character, for internal nodes,
		// it is the sum of the lower node's count
		unsigned int count;

		// indicates that this not a character node
		bool internal_node;

		// parent node of this node
		node* parent;

		// left child node
		node* left;

		// right child node
		node* right;
	};

	// human bit encoding
	struct code
	{
		int length;
		unsigned int value;
	};

	// cursor indicating the current point in a bit stream
	struct bit_cursor
	{
		unsigned char* bytePtr;
		unsigned char curBit;
	};

	// cursor indicating the current point in a read-only bit stream
	struct bit_cursor_read
	{
		const unsigned char* bytePtr;
		unsigned char curBit;
	};

	// pool of available nodes to be used as internal nodes of the huffman tree
	node internal_nodes[512];

	// array of leaf nodes, terminated by a sentinal with a maximum count value
	node sorted_leaf_nodes[256 + 1];

	// head of the huffman tree
	node* head = nullptr;

	// the character count of the input file, indexed by character
	unsigned int char_frequency[256] = {};

	// the maximum frequency in char_frequency
	unsigned int max_freq = 0;

	// the generated huffman codes, indexed by character
	code codes[256];

	// constant used to indicate the sentinal node at the end fo the leaf node list
	const unsigned int LAST_NODE = UINT_MAX;

	// sorts an input array by a specific decimal place
	void sort_radix(int r, unsigned char* in_array, unsigned char* out_array)
	{
		int div = pow(10, r);

		// count of each symbol (0-9)
		int symbol_count[10] = {};

		for (int i = 0; i < 256; i++)
		{
			symbol_count[(char_frequency[in_array[i]] / div) % 10]++;
		}

		// generate index offset table
		int offsets[10] = {};
		for (int i = 1; i < 10; i++)
		{
			offsets[i] = offsets[i - 1] + symbol_count[i - 1];
		}

		// sort array
		for (int i = 0; i < 256; i++)
		{
			int symbol = (char_frequency[in_array[i]] / div) % 10;
			out_array[offsets[symbol]] = in_array[i];
			offsets[symbol]++;
		}
	}

	// sort characters by their frequency and generate the leaf nodes for the huffman tree
	void sort_chars()
	{
		// copy char_frequency into sorted_chars
		int max_radix = ceil(log10(max_freq));

		unsigned char sorting_arrays[2][256];
		int input_array = 0;

		for (int i = 0; i < 256; i++)
		{
			sorting_arrays[input_array][i] = i;
		}

		for (int i = 0; i < max_radix; i++)
		{
			sort_radix(i, sorting_arrays[input_array], sorting_arrays[!input_array]);
			input_array = !input_array;
		}

		// output in descending order
		int leaf_idx = 0;
		for (int i = 0; i < 256; i++)
		{
			unsigned char character = sorting_arrays[input_array][i];
			if (char_frequency[character] > 0)
			{
				sorted_leaf_nodes[leaf_idx++] = { character, char_frequency[character], false, nullptr, nullptr, nullptr };
			}
		}

		// put a maximum-valued sentinal at the end of the list
		sorted_leaf_nodes[leaf_idx] = { '\0', LAST_NODE, false, nullptr, nullptr, nullptr };
	}

	// generates the huffman codes from the huffman tree
	void gen_codes(node* head, code cur_code)
	{
		if (head->internal_node)
		{
			cur_code.length++;
			cur_code.value <<= 1;
			gen_codes(head->left, cur_code);
			cur_code.value |= 1;
			gen_codes(head->right, cur_code);
		}
		else
		{
			codes[head->character] = cur_code;
		}
	}

	// writes bits with the given value (and "bits" being the number of bits) to the location
	// at the bit cursor
	void bitWrite(int bits, int value, bit_cursor& cursor)
	{
		while (bits > 0)
		{
			int bitsLeftInByte = 8 - cursor.curBit;
			if (bits < bitsLeftInByte)
			{
				*cursor.bytePtr |= value << (bitsLeftInByte - bits);
				cursor.curBit += bits;
				if (cursor.curBit >= 8)
				{
					cursor.curBit = 0;
					cursor.bytePtr++;
				}
				bits = 0;
			}
			else
			{
				*cursor.bytePtr |= (value >> (bits - bitsLeftInByte));
				value &= ((1 << (bits - bitsLeftInByte)) - 1);
				cursor.curBit = 0;
				cursor.bytePtr++;
				bits -= bitsLeftInByte;
			}
		}
	}

	// reads bits with the given number of bits (given by the parameter "bits") to the location
	// at the bit cursor
	unsigned int bitRead(int bits, bit_cursor_read& cursor)
	{
		auto get_bits = [&cursor](int bits) mutable -> unsigned char
		{
			unsigned char mask = ((1 << bits) - 1);
			return (*cursor.bytePtr >> (8 - cursor.curBit - bits)) & mask;
		};

		unsigned int data = 0;
		while (bits > 0)
		{
			int bitsLeftInByte = 8 - cursor.curBit;
			if (bitsLeftInByte > bits)
			{
				unsigned int dataRead = get_bits(bits);
				cursor.curBit += bits;
				data <<= bits;
				data |= dataRead;
				bits = 0;
			}
			else
			{
				unsigned int dataRead = get_bits(bitsLeftInByte);
				bits -= bitsLeftInByte;
				cursor.bytePtr++;
				cursor.curBit = 0;
				data <<= bitsLeftInByte;
				data |= dataRead;
			}
		}
		return data;
	}

public:

	// counts the frequency of each character type in the input buffer
	void count_chars(const unsigned char *bufin,
		unsigned int bufinlen)
	{
		max_freq = 0;
		for (unsigned int i = 0; i < bufinlen; i++)
		{
			char_frequency[bufin[i]]++;
			if (char_frequency[bufin[i]] > max_freq)
			{
				max_freq = char_frequency[bufin[i]];
			}
		}
	}

	// generate the huffman tree
	void huff_tree()
	{
		// sort characters by frequency and generate leaf nodes
		sort_chars();

		// priority queue of internal nodes sorted by count
		auto cmp = [](node* left, node* right) {
			return (left->count) > (right->count);
		};
		std::priority_queue<node*, std::vector<node*>, decltype(cmp)> node_q(cmp);

		// indices for node buffers
		int leaf_idx = 0;
		int internal_idx = 0;

		// finds the node with the lower count between leaf nodes and the internal nodes
		auto get_min_node = [this, &node_q, &leaf_idx]() mutable -> node*
		{
			if (node_q.empty() || sorted_leaf_nodes[leaf_idx].count < node_q.top()->count)
			{
				return &sorted_leaf_nodes[leaf_idx++];
			}
			else
			{
				node* ret = node_q.top();
				node_q.pop();
				return ret;
			}
		};
		
		// builds a new internal node with the given left and right nodes
		auto build_node = [this, &internal_idx](node* left, node* right) mutable -> node*
		{
			node* new_node = &internal_nodes[internal_idx++];
			new_node->internal_node = true;
			new_node->left = left;
			new_node->right = right;
			new_node->count = left->count + right->count;
			new_node->parent = nullptr;

			left->parent = new_node;
			right->parent = new_node;
			return new_node;
		};

		head = nullptr;
		while (head == nullptr)
		{
			node *left, *right;

			left = get_min_node();
			if (sorted_leaf_nodes[leaf_idx].count == LAST_NODE && node_q.empty())
			{
				// we've processed every node
				head = left;
			}
			else
			{
				right = get_min_node();
				node_q.push(build_node(left, right));
			}
		}
	}

	// generates the huffman codes from the huffman tree
	void gen_codes()
	{
		gen_codes(head, { 0, 0 });
	}

	// calculates the size of the output buffer
	unsigned int calc_out_buff_size()
	{
		unsigned long bit_count = 0;
		for (int i = 0; i < 256; i++)
		{
			bit_count += char_frequency[i] * codes[i].length;
		}
		return ((bit_count + 7) / 8) + sizeof(char_frequency);
	}

	// writes the character frequency table to the output file
	void write_char_freq(unsigned char*& pbufout)
	{
		unsigned int* intPtr = (unsigned int *)pbufout;
		for (int i = 0; i < 256; i++)
		{
			*intPtr = char_frequency[i];
			intPtr++;
		}
		pbufout = (unsigned char *)intPtr;
	}

	// writes the huffman encoding to the output file
	void encode(const unsigned char *bufin,
		unsigned int bufinlen,
		unsigned char *pbufout)
	{
		bit_cursor outCursor = { pbufout, 0 };
		for (unsigned int i = 0; i < bufinlen; i++)
		{
			bitWrite(codes[bufin[i]].length, codes[bufin[i]].value, outCursor);
		}
	}

	// read the character frequency table from the input file
    void read_char_freq(const unsigned char*& bufin, unsigned int& bufinlen, 
		unsigned int& totalSize, unsigned int& uniqueChars)
	{
		uniqueChars = 0;
		max_freq = 0;
		unsigned int* intPtr = (unsigned int *)bufin;
		for (int i = 0; i < 256; i++)
		{
			char_frequency[i] = *intPtr;
			totalSize += char_frequency[i];
			if (char_frequency[i] > max_freq)
			{
				max_freq = char_frequency[i];
			}
			if (char_frequency[i] > 0)
			{
				uniqueChars++;
			}
			intPtr++;
		}
		bufinlen -= sizeof(char_frequency);
		bufin = (unsigned char *)intPtr;
	}

	// decodes the input file until all the characters in the character
	// frequency table are accounted for
	void decode(const unsigned char *bufin,
		unsigned int bufinlen,
		unsigned char *pbufout,
		unsigned int uniqueChars)
	{
		bit_cursor_read inCursor = { bufin, 0 };
		node* currentNode = head;
		while (uniqueChars > 0)
		{
			unsigned char curBit = bitRead(1, inCursor);

			// NOTE: this check is needed for the case where there is only one type of character so the
			// head of the tree would be a leaf node
			if (currentNode->internal_node)
			{
				currentNode = curBit == 0 ? currentNode->left : currentNode->right;
			}

			if (!currentNode->internal_node)
			{
				*pbufout = currentNode->character;
				pbufout++;
				char_frequency[currentNode->character]--;
				if (char_frequency[currentNode->character] == 0)
				{
					uniqueChars--;
				}
				currentNode = head;
			}
		}
	}
};

int huffman_encode(const unsigned char *bufin,
	unsigned int bufinlen,
	unsigned char **pbufout,
	unsigned int *pbufoutlen)
{
	huffman_codec codec;

	// build frequency table
	codec.count_chars(bufin, bufinlen);

	// build huffman tree
	codec.huff_tree();

	// generate huffman codes
	codec.gen_codes();

	// allocate output buffer
	*pbufoutlen = codec.calc_out_buff_size();
	*pbufout = (unsigned char *)malloc(*pbufoutlen);
	unsigned char* curPtr = *pbufout;
	memset(curPtr, 0, *pbufoutlen);

	// output character frequency table
	codec.write_char_freq(curPtr);

	// encode input buffer
	codec.encode(bufin, bufinlen, curPtr);
	return 0;
}

int huffman_decode(const unsigned char *bufin,
	unsigned int bufinlen,
	unsigned char **pbufout,
	unsigned int *pbufoutlen)
{
	unsigned int uniqueChars;
	huffman_codec codec;

	// read character frequency table
	codec.read_char_freq(bufin, bufinlen, *pbufoutlen, uniqueChars);

	// build huffman tree
	codec.huff_tree();

	// generate output buffer (leave space for null terminator needed for output stream to file)
	*pbufout = (unsigned char *)malloc(*pbufoutlen + 1);

	// decode input buffer
	codec.decode(bufin, bufinlen, *pbufout, uniqueChars);

	// append null terminator
	(*pbufout)[*pbufoutlen] = '\0';
	*pbufoutlen = *pbufoutlen+1;

	return 0;
}